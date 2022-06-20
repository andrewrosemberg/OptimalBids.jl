using OptimalBids
using OptimalBids.PowerModelsMarkets
using Gurobi # Market Clearing Solver
using JuMP: optimizer_with_attributes

using Nonconvex
using NonconvexIpopt # Nonconvex.@load Ipopt
using NonconvexBayesian # Nonconvex.@load BayesOpt
using Zygote

using AbstractGPs
using KernelFunctions

using Plots # For some evaluation plots at the end
using Plots.PlotMeasures
using Downloads # To download Test Cases

# # Case Definition 

# Read market data from IEEE 118 bus case
case_name = "pglib_opf_case118_ieee.m"
DATA_DIR = mktempdir()
case_file_path = joinpath(DATA_DIR, case_name)
Downloads.download(
    "https://raw.githubusercontent.com/power-grid-lib/pglib-opf/01681386d084d8bd03b429abcd1ee6966f68b9a3/" *
    case_name,
    case_file_path,
)
network_data = PowerModels.parse_file(case_file_path)

# Measure maximum load
max_load =
    sum(load["pd"] for load in values(network_data["load"])) * network_data["baseMVA"]

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

# # OptimalBids API
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
offer_weights = offer_weights / sum(offer_weights)

# However, the decision maker is allowed to increase all bids evenly:
min_total_volume = 0.0
max_total_volume = 155.0
range_mul_factor = min_total_volume:1.0:max_total_volume
bid_range = [offer_weights .* [i] for i in range_mul_factor]
p_curve = profit_curve!(market, bid_range)
maximum_pq_curve, argmax_pq_curve = findmax(p_curve)

# Mesure time to change bids and solve opf
num_opfs = 10
opf_time = @elapsed [profit_for_bid!(market, offer_weights) for _ in 1:num_opfs]
opf_time /= num_opfs

# Let's plot and see how the range profit evaluatiuon went:
plt_range = plot(
    collect(range_mul_factor) * 100 / max_load,
    p_curve;
    label="Range Evaluation",
    ylabel="Profit (USD)",
    xlabel="Bid Volume (% Market Share)",
    legend=:outertopright,
    left_margin=10mm,
    bottom_margin=10mm,
    size=(900, 600),
    color=:black,
    width=2.0,
);

# # Nonconvex API

# Nonconvex needs a minimization objective function that only receives the decision vector.
function profit_function(total_volume)
    return -profit_for_bid!(market, offer_weights .* first(total_volume))
end

mutable struct StorageCallbackProfit <: Function
    start_time::Float64
    fcalls::Int
    visited_objective::Array{Float64}
    visited_volumes::Array{Float64}
    visited_times::Array{Float64}
    projection::Array{Float64}
end
function StorageCallbackProfit(maxiter, start_time; projection=offer_weights)
    return StorageCallbackProfit(
        start_time,
        0,
        Array{Float64}(undef, 5 * maxiter),
        Array{Float64}(undef, 5 * maxiter),
        Array{Float64}(undef, 5 * maxiter),
        projection,
    )
end

function (callback::StorageCallbackProfit)(total_volume)
    Zygote.@ignore callback.fcalls += 1
    Zygote.@ignore callback.visited_volumes[callback.fcalls] = sum(
        callback.projection .* total_volume
    )
    Zygote.@ignore callback.visited_times[callback.fcalls] = time() - callback.start_time
    obj = profit_for_bid!(market, callback.projection .* total_volume)
    Zygote.@ignore callback.visited_objective[callback.fcalls] = obj
    return -obj
end

# # BayesOpt (0-Order)

# Max Number of Iterations for the solution method (proxy to a time limit at bidding time).
# ps.: Currently, no option for limiting fcalls.
maxiter = 10

storage_profit_function = StorageCallbackProfit(maxiter, time())

# Build Nonconvex optimization model:
model = Nonconvex.Model()
set_objective!(model, storage_profit_function; flags=[:expensive])
addvar!(model, [min_total_volume], [max_total_volume])
add_ineq_constraint!(model, x -> sum(x) - max_total_volume)

# Solution Method: Bayesian Optimization
alg = BayesOptAlg(IpoptAlg())
options = BayesOptOptions(;
    sub_options=IpoptOptions(; max_iter=20, print_level=0),
    # ninit=Int(floor(maxiter / 5)),
    maxiter=maxiter,
    ftol=1e-4,
    ctol=1e-5,
    initialize=true,
    postoptimize=false,
    kernel=RationalKernel(; α=2.27e8) ∘ ScaleTransform(0.01),
    noise=0.001,
    std_multiple=8.67e4,
    fit_prior=false, # not working with custom priors
)

# Optimize model:
r_bayes = optimize(model, alg, [max_total_volume / 2]; options=options)

best_solution = r_bayes.minimizer
best_profit = -r_bayes.minimum

x = storage_profit_function.visited_volumes[1:(storage_profit_function.fcalls)]
y = storage_profit_function.visited_objective[1:(storage_profit_function.fcalls)]

plt_gif = Animation()
for iter in 1:(storage_profit_function.fcalls)
    # Final GP
    f = GP(8.67e8 * RationalKernel(; α=2.27e8) ∘ ScaleTransform(0.01))

    # Finite projection of `f` at inputs `x`.
    # Added Gaussian noise with variance 0.001.
    fx = f(x[1:iter], 0.001)

    # Exact posterior given `y`. This is another GP.
    p_fx = posterior(fx, y[1:iter])

    # Plot posterior.
    plt_range_aux = plot(
        collect(range_mul_factor) * 100 / max_load,
        p_curve;
        label="Range Evaluation",
        ylabel="Profit (USD)",
        xlabel="Bid Volume (% Market Share)",
        title="BayesOpt Analysis",
        legend=:outertopright,
        left_margin=10mm,
        bottom_margin=10mm,
        size=(900, 600),
        color=:black,
        width=2.0,
    )
    plot!(
        plt_range_aux,
        collect(range_mul_factor) * 100 / max_load,
        p_fx(collect(range_mul_factor));
        label="Posterior - " * string(RationalKernel),
    )
    scatter!(plt_range_aux, x[1:iter] * 100 / max_load, y[1:iter]; label="Visited")
    frame(plt_gif, plt_range_aux)
end

gif(plt_gif, joinpath(pwd(), "examples/PowerModelsMarkets/GaussianProcesses/GP.gif"); fps=3)
