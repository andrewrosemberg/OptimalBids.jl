# Profit Maximization Example
## using Nonconvex.jl
This a example on how to use [Nonconvex.jl](https://github.com/JuliaNonconvex/Nonconvex.jl) to maximize the profit of a company operating in a defined market (through the `OptimalBids` API).

 - For this example, we will use a market of type `PowerModelsMarkets` which represents an energy spot market.
 - Case instance is IEEE 118 Bus Case
 - To ilustrate a situation where the company defines the distribution of volume across it's assets in a predefined way (and reduce the decision dimentionality), the only controlable variable to maximize the company profit is a multiplicative factor of all offers.  

### Include Needed Packages

```@example Nonconvex
using OptimalBids
using OptimalBids.PowerModelsMarkets
using Gurobi # Market Clearing Solver
using JuMP: optimizer_with_attributes

using Nonconvex
using NonconvexIpopt # Nonconvex.@load Ipopt
using NonconvexBayesian # Nonconvex.@load BayesOpt

using AbstractGPs
using KernelFunctions

using Plots # For some evaluation plots at the end
using Plots.PlotMeasures
using Downloads # To download Test Cases

```

### Case Definition

```@example Nonconvex

# Read market data from IEEE 118 bus case
case_name = "pglib_opf_case118_ieee.m"
DATA_DIR = mktempdir()
case_file_path = joinpath(DATA_DIR, case_name)
Downloads.download("https://raw.githubusercontent.com/power-grid-lib/pglib-opf/01681386d084d8bd03b429abcd1ee6966f68b9a3/" * case_name, case_file_path)
network_data = PowerModels.parse_file(case_file_path)

# Pretend we are a company constructing a new set of generators in the grid.
# Choose a percentage of the total number of buses to install the new generators:
percentage_buses = 0.09

# We need the keys PowerModels uses to reference the appropriate buses in it's network data dictionary.
# First, find out all available keys:
bus_indexes = collect(keys(network_data["bus"]))
# Then, calculate number of buses that consitute the chose percent (`percentage_buses`):
num_buses = length(bus_indexes)
num_strategic_buses = ceil(Int, percentage_buses * num_buses)
# To avoid any biases let's grab some generators in the middle:
bus_indexes = bus_indexes[21:(21 + num_strategic_buses - 1)]
# Finally, add new generators to the network grid data and collect their reference keys.
generator_indexes = [
    add_generator(network_data, parse(Int, bus_idx)) for bus_idx in bus_indexes
]

```

### OptimalBids API

```@example Nonconvex

# Define market
market = build_market(
    PowerModelsMarket,
    network_data,
    generator_indexes,
    optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0),
)

# Relative distribution of offers are sometimes predefined and cannot be changed at bidding time.
using Random
rng = MersenneTwister(666)
offer_weights = rand(rng, num_strategic_buses)
offer_weights = offer_weights/ sum(offer_weights)

# However, the decision maker is allowed to increase all bids evenly:
min_total_volume = 0.0
max_total_volume = 150.0
range_mul_factor = min_total_volume:1.0:max_total_volume
bid_range = [offer_weights .* [i] for i in range_mul_factor]
p_curve = profit_curve!(market, bid_range)
maximum_pq_curve, argmax_pq_curve = findmax(p_curve)

# Let's plot and see how the range profit evaluatiuon went:
plt_range = plot(collect(range_mul_factor), p_curve,
    label="Range Evaluation",
    ylabel="Profit (\$)",
    xlabel="Total Volume (MWh)",
    legend=:outertopright,
    left_margin=10mm,
    bottom_margin=10mm,
    size=(900, 600)
);
plt_comp = deepcopy(plt_range);
```

### Nonconvex API

```@example Nonconvex

# Nonconvex needs a minimization objective function that only receives the decision vector.
function profit_function(total_volume)
    return - profit_for_bid!(market, offer_weights .* first(total_volume))
end

mutable struct StorageCallbackProfit <: Function
    start_time::Float64
    fcalls::Int
    visited_objective::Array{Float64}
    visited_volumes::Array{Float64}
    visited_times::Array{Float64}
end
StorageCallbackProfit(maxiter, start_time) = StorageCallbackProfit(start_time, 0, 
    Array{Float64}(undef, 5 * maxiter),
    Array{Float64}(undef, 5 * maxiter),
    Array{Float64}(undef, 5 * maxiter),
)

function (callback::StorageCallbackProfit)(total_volume)
    Zygote.@ignore callback.fcalls += 1
    Zygote.@ignore callback.visited_volumes[callback.fcalls] = first(total_volume)
    Zygote.@ignore callback.visited_times[callback.fcalls] = time() - callback.start_time
    obj = profit_for_bid!(market, offer_weights .* first(total_volume))
    Zygote.@ignore callback.visited_objective[callback.fcalls] = obj
    return - obj
end

# Find time to change bids and solve opf
num_opfs = 10
opf_time = @elapsed [profit_for_bid!(market, offer_weights) for _ = 1:num_opfs]
opf_time /= num_opfs

```

### BayesOpt (0-Order)

```@example Nonconvex

# Max Number of Iterations for the solution method (proxy to a time limit at bidding time).
# ps.: Currently, no option for limiting fcalls.
maxiter = 10

storage_profit_function = StorageCallbackProfit(maxiter, time())

# Build Nonconvex optimization model:
model = Nonconvex.Model()
set_objective!(model, storage_profit_function, flags = [:expensive])
addvar!(model, [min_total_volume], [max_total_volume])
add_ineq_constraint!(model, x -> -1) # errors when no inequality is added!

# Solution Method: Bayesian Optimization
alg = BayesOptAlg(IpoptAlg())
options = BayesOptOptions(
    sub_options = IpoptOptions(max_iter = 20, print_level = 0),
    # ninit=Int(floor(maxiter / 5)),
    maxiter = maxiter, ftol = 1e-4, ctol = 1e-5, initialize=true, postoptimize=false,
    kernel= RationalKernel(α=2.27e8) ∘ ScaleTransform(0.01),
    noise=0.001,
    std_multiple=8.67e4,
    fit_prior=false # not working with custom priors
)

# Optimize model:
r_bayes = optimize(model, alg, [max_total_volume / 2]; options = options)

best_solution = r_bayes.minimizer
best_profit = -r_bayes.minimum

scatter!(plt_comp, [best_solution], [best_profit],
    label="BayesOpt - OPF Calls:$(storage_profit_function.fcalls)",
)

plt_surrogate = deepcopy(plt_range);

up_surrugate = -getproperty.(r_bayes.surrogates[1].(range_mul_factor), :lo)
lb_surrugate = -getproperty.(r_bayes.surrogates[1].(range_mul_factor), :hi)
std_surrugate = (up_surrugate .- lb_surrugate) / 2
med_surrugate = lb_surrugate + std_surrugate

storage_profit_function.visited_times = storage_profit_function.visited_times ./ opf_time

plot!(plt_surrogate, range_mul_factor, med_surrugate,
    ribbon=std_surrugate,
    title="BayesOpt Analysis",
    label="Surrogate Function",
);
scatter!(plt_surrogate, storage_profit_function.visited_volumes[1:storage_profit_function.fcalls], 
    storage_profit_function.visited_objective[1:storage_profit_function.fcalls]; label="Visited"
)
```

```@example Nonconvex

plot(storage_profit_function.visited_times[1:storage_profit_function.fcalls], 
(maximum_pq_curve .- accumulate(max, storage_profit_function.visited_objective[1:storage_profit_function.fcalls])) ./ maximum_pq_curve,
    xlabel="Time (x OPF)",
    ylabel="Optimality Gap (%)",
    ylim=(0.0,1.0),
    legend=false
)
```
### NLopt - BOBYQA (0-Order)

```@example Nonconvex

using NonconvexNLopt

maxeval = 50
storage_profit_function = StorageCallbackProfit(maxeval, time())

# Build Nonconvex optimization model:
model = Nonconvex.Model()
set_objective!(model, storage_profit_function)
addvar!(model, [min_total_volume], [max_total_volume])

# Solution Method: Sequential Least-Squares Quadratic Programming
method = :LN_BOBYQA
alg = NLoptAlg(:LN_BOBYQA)
options = NLoptOptions(maxeval=maxeval)

# Optimize model
callback_storage = StorageCallback(maxeval, time())
r = optimize(model, alg, [max_total_volume / 2]; options = options)

best_solution = r.minimizer
best_profit = -r.minimum

scatter!(plt_comp, [best_solution], [best_profit],
    label="NLOpt-$(method) - OPF Calls:$(storage_profit_function.fcalls)",
);


plt_surrogate = deepcopy(plt_range);

storage_profit_function.visited_times = storage_profit_function.visited_times ./ opf_time

scatter!(plt_surrogate, storage_profit_function.visited_volumes[1:storage_profit_function.fcalls], 
    storage_profit_function.visited_objective[1:storage_profit_function.fcalls]; label="Visited"
)

plot(storage_profit_function.visited_times[1:storage_profit_function.fcalls], 
(maximum_pq_curve .- accumulate(max, storage_profit_function.visited_objective[1:storage_profit_function.fcalls])) ./ maximum_pq_curve,
    xlabel="Time (x OPF)",
    ylabel="Optimality Gap (%)",
    ylim=(0.0,1.0),
    legend=false
)
```

### NonconvexMultistart - GPSampler (0-Order)

```@example Nonconvex

using NonconvexMultistart

# Build Nonconvex optimization model:
model = Nonconvex.Model()
set_objective!(model, profit_function)
addvar!(model, [min_total_volume], [max_total_volume])

# Solution Method: Hyperopt
global fcalls = 0
maxiter = 4
method = :Hyperopt
alg = HyperoptAlg(IpoptAlg())
options = HyperoptOptions(
    sub_options = IpoptOptions(max_iter = maxiter), sampler = GPSampler(),
    iters = 2,
    keep_all=true
)

# Optimize model
r_hyp = optimize(model, alg, [min_total_volume], options = options)

best_solution = r_hyp.minimizer
best_profit = -r_hyp.minimum

scatter!(plt_comp, best_solution, [best_profit],
    label="$(method) Offer - OPF Calls:$(fcalls)",
);

plt_visited = deepcopy(plt_range)

scatter!(plt_visited, [i.minimizer[1] for i in r_hyp.results], [- i.minimum for i in r_hyp.results],
    label="$(method) Offer - OPF Calls:$(fcalls)",
    title="Hyperopt Visited Offers",
)
```

### Profit Comparison NLP 0-Order Strategies

```@example Nonconvex

plot(plt_comp, margin=5Plots.mm,
    title="Profit Comparison NLP Strategies",
)
```

### NLopt - CCSAQ (1-Order)

```@example Nonconvex

using NonconvexNLopt

maxeval = 4
global fcalls = 0

# Build Nonconvex optimization model:
model = Nonconvex.Model()
set_objective!(model, profit_function)
addvar!(model, [min_total_volume], [max_total_volume])

# Solution Method: Sequential Least-Squares Quadratic Programming
method = :LD_CCSAQ
alg = NLoptAlg(method)
options = NLoptOptions(maxeval=maxeval)

# Optimize model
r = optimize(model, alg, [min_total_volume], options = options)

best_solution = r.minimizer
best_profit = -r.minimum

scatter!(plt_comp, [best_solution], [best_profit],
    label="NLOpt-$(method) - OPF Calls:$(fcalls)",
)
```