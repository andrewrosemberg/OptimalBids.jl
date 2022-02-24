## Profit Maximization Example: using Nonconvex.jl
This a example on how to use [Nonconvex.jl](https://github.com/JuliaNonconvex/Nonconvex.jl) to maximize the profit of a company operating in a defined market (through the `OptimalBids` API).

For this example, we will use a market of type `PowerModelsMarkets`. 

```@example Nonconvex
using OptimalBids
using OptimalBids.PowerModelsMarkets
using Clp # Market Clearing Solver

using Nonconvex
using NonconvexIpopt # Nonconvex.@load Ipopt
using NonconvexBayesian # Nonconvex.@load BayesOpt

using Plots # For some evaluation plots at the end

#=
CASE DEFINITION
=#

# Read market data from IEEE 118 bus case
case_name = "case118.m"
DATA_DIR = joinpath(dirname(dirname(dirname(@__DIR__))), "test/data")
case_file_path = joinpath(DATA_DIR, case_name)
network_data = PowerModels.parse_file(case_file_path)

# Let's make the case a bit more interesting, by adding some randomness to existing generators costs and available load.
using Random
Random.seed!(654654)
for gen in values(network_data["gen"])
    gen["cost"][end-1] += rand(-2000.0:2500.0)
end
load_mul_factor = 6.8
for load in collect(values(network_data["load"]))[1:20:end]
    load["pd"] *= load_mul_factor * rand(0.01:0.01:3)
end

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

#=
OptimalBids API
=#

# Define market
market = build_market(
    PowerModelsMarket,
    network_data,
    generator_indexes,
    Clp.Optimizer,
)

# Relative distribution of offers are sometimes predefined and cannot be changed bidding time.
offer_weights = rand(num_strategic_buses)
offer_weights = offer_weights/ sum(offer_weights)

# However, the decision maker is allowed to increase all bids evenly:
min_total_volume = 0.0
max_total_volume = 65.0
range_mul_factor = min_total_volume:0.1:max_total_volume
bid_range = [offer_weights .* [i] for i in range_mul_factor]
p_curve = profit_curve!(market, bid_range)

# Let's plot and see how the range profit evaluatiuon went:
plt_range = plot(collect(range_mul_factor), p_curve,
    title="Case $case_name",
    label="Range Evaluation - Random Offers",
    ylabel="Profit (\$)",
    xlabel="Multiplicative Factor",
    legend=:outertopright,
)

#=
Nonconvex API
=#

## BayesOpt

# Nonconvex needs a minimization objective function that only receives the decision vector.
function profit_function(total_volume)
    return - profit_for_bid!(market, offer_weights .* total_volume[1])
end

# Max Number of Iterations for the solution method (proxy to a time limit at bidding time).
# ps.: Currently, no option for limiting fcalls.
maxiter = 10

# Build Nonconvex optimization model:
model = Model()
set_objective!(model, profit_function, flags = [:expensive])
addvar!(model, [min_total_volume], [max_total_volume])
add_ineq_constraint!(model, x -> -1) # errors when no inequality is added! 

# Solution Method: Bayesian Optimization
alg = BayesOptAlg(IpoptAlg())
options = BayesOptOptions(
    sub_options = IpoptOptions(max_iter = maxiter),
    maxiter = maxiter, ftol = 1e-4, ctol = 1e-5,
)

# Optimize model:
r_bayes = optimize(model, alg, [min_total_volume], options = options)

best_solution = r_bayes.minimizer
best_profit = -r_bayes.minimum
r_bayes.niters # number of iterations of the 

scatter!(plt_range, [best_solution; r_bayes.sub_result.minimizer], [best_profit; -r_bayes.sub_result.minimum],
    label="BayesOpt Offer - OPF Calls:$(r_bayes.sub_result.fcalls)",
    size=(1000, 1000)
)

plt_surrogate = deepcopy(plt_range)

plot!(plt_surrogate, range_mul_factor, -getproperty.(r_bayes.surrogates[1].(range_mul_factor), :lo),
    label="BayesOpt - Surrogate Function",
)

## NLOpt
using NonconvexNLopt

include(joinpath(dirname(dirname(dirname(@__DIR__))), "test/fix_nlopt.jl")) # Issue: https://github.com/andrewrosemberg/OptimalBids.jl/issues/8

maxeval = 10

# Build Nonconvex optimization model:
model = Model()
set_objective!(model, profit_function)
addvar!(model, [min_total_volume], [max_total_volume])

# Solution Method: Bayesian Optimization
method = :LN_BOBYQA
alg = NLoptAlg(:LN_BOBYQA)
options = NLoptOptions(maxeval=maxeval)

# Optimize model: Sequential Least-Squares Quadratic Programming
r = optimize(model, alg, [min_total_volume], options = options)

best_solution = r.minimizer
best_profit = -r.minimum
r.fcalls # number of function calls

scatter!(plt_range, [best_solution], [best_profit],
    label="NLOpt-$(method) Offer - OPF Calls:$(r.fcalls)",
    size=(1000, 1000)
)

plot(plt_range, plt_surrogate, size=(1400, 1000))

```

![](https://github.com/andrewrosemberg/PortfolioOpt/blob/master/docs/src/assets/bayesopt_profit.png?raw=true)