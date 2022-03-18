# Examples
## Profit Maximization Example: using Optim.jl
This a example on how to use [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl/) to maximize the profit of a company operating in a defined market (through the `OptimalBids` API).

 - For this example, we will use a market of type `PowerModelsMarkets` which represents an energy spot market.
 - Case instance is IEEE 118 Bus Case
 - To ilustrate a situation where the company defines the distribution of volume across it's assets in a predefined way (and reduce the decision dimentionality), the only controlable variable to maximize the company profit is a multiplicative factor of all offers.  

### Include Needed Packages

```@example Optim
using OptimalBids
using OptimalBids.PowerModelsMarkets
using Clp # Market Clearing Solver
using JuMP: optimizer_with_attributes

using Optim

using Plots # For some evaluation plots at the end
using Downloads # To download Test Cases

```

### Case Definition

```@example Optim

# Read market data from IEEE 118 bus case
case_name = "pglib_opf_case118_ieee.m"
DATA_DIR = mktempdir()
case_file_path = joinpath(DATA_DIR, case_name)
Downloads.download(joinpath("https://raw.githubusercontent.com/power-grid-lib/pglib-opf/01681386d084d8bd03b429abcd1ee6966f68b9a3", case_name), case_file_path)
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

```@example Optim

# Define market
market = build_market(
    PowerModelsMarket,
    network_data,
    generator_indexes,
    optimizer_with_attributes(Clp.Optimizer, "LogLevel" => 0),
)

# Relative distribution of offers are sometimes predefined and cannot be changed at bidding time.
using Random
rng = MersenneTwister(0)
offer_weights = rand(rng, num_strategic_buses)
offer_weights = offer_weights/ sum(offer_weights)

# However, the decision maker is allowed to increase all bids evenly:
min_total_volume = 0.0
max_total_volume = 65.0
range_mul_factor = min_total_volume:0.1:max_total_volume
bid_range = [offer_weights .* [i] for i in range_mul_factor]
p_curve = profit_curve!(market, bid_range)

# Let's plot and see how the range profit evaluatiuon went:
plt_range = plot(collect(range_mul_factor), p_curve,
    label="Range Evaluation",
    ylabel="Profit (\$)",
    xlabel="Multiplicative Factor",
    legend=:outertopright,
);
plt_comp = deepcopy(plt_range);
```

### Profit Function

```@example Optim

# Optim needs a minimization objective function that only receives the decision vector.
function profit_function(total_volume)
    global fcalls += 1
    return - profit_for_bid!(market, offer_weights .* total_volume[1])
end

```

### NelderMead (0-Order)

```@example Optim

# Max Number of Iterations for the solution method (proxy to a time limit at bidding time).
maxiter = 10
global fcalls = 0

# Solve Optim optimization model:
lower = [min_total_volume]
upper = [max_total_volume]
initial_x = [min_total_volume] .+ 0.001 # Initial value cannot be on boundary value

# Limiting fcalls.
function advanced_fcall_control(x)
    println(" * Fcalls so far:     ", fcalls)
    if fcalls >= maxiter
        return true
    end
    println()
    false
end

inner_optimizer = NelderMead()
results = Optim.optimize(profit_function, lower, upper, initial_x, Fminbox(inner_optimizer), Optim.Options(f_calls_limit=maxiter, callback = advanced_fcall_control, store_trace=true))

best_solution = results.minimizer
best_profit = -results.minimum

scatter!(plt_comp, [best_solution], [best_profit],
    label="NelderMead - OPF Calls:$(fcalls)",
)
```

### LBFGS (1-Order)

```@example Optim

# Profit and gradient
using LinearAlgebra

function fg!(F,G,total_volume)
    y, Δ = profit_and_gradient_for_bid!(market, offer_weights .* total_volume[1])
    global fcalls += 1

    isnothing(G) || copyto!(G, dot(Δ, offer_weights))
    isnothing(F) || return y

    nothing
end

# Max Number of Iterations for the solution method (proxy to a time limit at bidding time).
maxiter = 10
global fcalls = 0

inner_optimizer = LBFGS()
results = Optim.optimize(Optim.only_fg!(fg!), lower, upper, initial_x, Fminbox(inner_optimizer), Optim.Options(f_calls_limit=maxiter, callback = advanced_fcall_control, store_trace=true))

best_solution = results.minimizer
best_profit = -results.minimum

scatter!(plt_comp, [best_solution], [best_profit],
    label="LBFGS - OPF Calls:$(fcalls)",
)
```

### Brent (0-Order)

```@example Optim

# Max Number of Iterations for the solution method (proxy to a time limit at bidding time).
maxiter = 10
global fcalls = 0

function profit_function(total_volume)
    global fcalls += 1
    return - profit_for_bid!(market, offer_weights .* total_volume)
end
results = Optim.optimize(profit_function, min_total_volume, max_total_volume, Brent(); callback = advanced_fcall_control, store_trace=true)

best_solution = results.minimizer
best_profit = -results.minimum

scatter!(plt_comp, [best_solution], [best_profit],
    label="Brent - OPF Calls:$(fcalls)",
)
```