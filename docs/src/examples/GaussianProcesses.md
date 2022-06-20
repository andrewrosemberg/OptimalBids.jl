# Example Fitting a Gaussian Process to the profit curve (using AbstractGPs.jl)

This a example on how to use [AbstractGPs.jl](https://github.com/JuliaGaussianProcesses/AbstractGPs.jl) to fit a Gaussian Process (GP) to the profit curve of a company operating in a defined market (through the `OptimalBids` API).

 - For this example, we will use a market of type `PowerModelsMarkets` which represents an energy spot market.
 - Case instance is IEEE 118 Bus Case
 - To ilustrate a situation where the company defines the distribution of volume across it's assets in a predefined way (and reduce the decision dimentionality), the only controlable variable to maximize the company profit is a multiplicative factor of all offers.
 - We will use [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl/) for the MLE of the tested kernels' hyperparameters. 

## Include Needed Packages


```@example AbstractGPs
using OptimalBids
using OptimalBids.PowerModelsMarkets
using Clp # Market Clearing Solver
using JuMP: optimizer_with_attributes

using AbstractGPs
using KernelFunctions

using Optim
using ParameterHandling
using Zygote
using ParameterHandling: flatten

using Plots # For some evaluation plots at the end
using Plots.PlotMeasures
using Downloads # To download Test Cases

```

## Case Definition

```@example AbstractGPs

# Read market data from IEEE 118 bus case
case_name = "pglib_opf_case118_ieee.m"
DATA_DIR = mktempdir()
case_file_path = joinpath(DATA_DIR, case_name)
Downloads.download("https://raw.githubusercontent.com/power-grid-lib/pglib-opf/01681386d084d8bd03b429abcd1ee6966f68b9a3/" * case_name, case_file_path)
network_data = PowerModels.parse_file(case_file_path)

# Measure maximum load
max_load = sum(load["pd"] for load in values(network_data["load"])) * network_data["baseMVA"]

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

## OptimalBids API

```@example AbstractGPs

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
max_total_volume = 655.0
range_mul_factor = min_total_volume:1.0:max_total_volume
bid_range = [offer_weights .* [i] for i in range_mul_factor]
p_curve = profit_curve!(market, bid_range)

```

## Optim API

```@example AbstractGPs

function bfgs(objective, initial_params; verbose=false)
   training_results = Optim.optimize(
        objective,
        θ -> only(Zygote.gradient(objective, θ)),
        initial_params,
        LBFGS(
            alphaguess = Optim.LineSearches.InitialStatic(scaled=true),
            linesearch = Optim.LineSearches.BackTracking(),
        ),
        Optim.Options(show_trace = verbose);
        inplace=false,
    )

   return (training_results.minimum, training_results.minimizer)
end

function multi_start_hyperopt(target, algorithm, num_iterations, initial_initial_params; update_initials = (θ_last, θ_best, obj_diff) -> rand(length(θ_last)), verbose=false, kwargs...)
    θ_initials = Array{Float64, 2}(undef, num_iterations, length(initial_initial_params));
    vistited_best_fs = Array{Float64}(undef, num_iterations);
    best_θ = fill(1/length(initial_initial_params), length(initial_initial_params))
    best_f = target(best_θ)
    last_θ = best_θ
    last_f = best_f
    for i = 1:num_iterations
        try
            (last_f, last_θ) = algorithm(target, update_initials(last_θ, best_θ, last_f-best_f); verbose=verbose, kwargs...)
        catch
            println("SHIT")
        end
        θ_initials[i,:] .= last_θ
        if last_f < best_f
            best_f = last_f
            best_θ = last_θ
        end
        vistited_best_fs[i] = best_f
    end
    return best_f, best_θ, θ_initials, vistited_best_fs
end

```

## AbstractGPs Interface

### Data

```@example AbstractGPs

# Randomly choose observed data
number_samples = 9
idx = [1; rand(2:size(range_mul_factor,1), number_samples)]
x = collect(range_mul_factor)[idx]
y = p_curve[idx]

plt_comp = plot(collect(range_mul_factor) * 100 / max_load, p_curve,
    label="Range Evaluation",
    ylabel="Profit (USD)",
    size=(900, 600), 
    title="Kernel Comparison", 
    xlabel="Bid Volume (% Market Share)",
    legend=:outertopright,
    color="black",
    width=3,
    left_margin=10mm,
    bottom_margin=10mm,
);
scatter!(plt_comp, x .* 100 ./ max_load, y; label="Data");

plt_comp
```

### Matern32Kernel

```@example AbstractGPs

## Define GP prior with Matern32Kernel
flat_initial_params, unflatten = flatten((;
    var_kernel = bounded(10.0, 1.0, 1e8),
    λ = bounded(0.005, 0.0, 0.01),
))

# Construct a function to unpack flattened parameters and pull out the raw values.
unpack = ParameterHandling.value ∘ unflatten
params = unpack(flat_initial_params)

function build_gp(params)
    return GP(params.var_kernel * Matern52Kernel() ∘ ScaleTransform(params.λ))
end

# Define MLE objective
function objective(params)
    f = build_gp(unpack(params))
    return -logpdf(f(x, 0.001), y)
end

# Optimise using Optim
(best_f, best_θ, θ_initials, vistited_best_fs) = multi_start_hyperopt(objective, bfgs, 100, flat_initial_params)

# Extracting the optimal values of the parameters
optimal_params = unpack(best_θ)

# Final GP
f = build_gp(optimal_params)

# Finite projection of `f` at inputs `x`.
# Added Gaussian noise with variance 0.001.
fx = f(x, 0.001)

# Exact posterior given `y`. This is another GP.
p_fx = posterior(fx, y)

# Plot posterior.
plot!(plt_comp, collect(range_mul_factor) * 100 / max_load, p_fx(collect(range_mul_factor))(collect(range_mul_factor)); label="Posterior - " * string(Matern32Kernel))
```


### RationalKernel

```@example AbstractGPs

## Define GP prior with RationalKernel
flat_initial_params, unflatten = flatten((;
    var_kernel = bounded(10.0, 6.0, 1e8),
    λ = bounded(0.005, 0.0, 0.01),
    α = positive(0.05),
))

# Construct a function to unpack flattened parameters and pull out the raw values.
unpack = ParameterHandling.value ∘ unflatten
params = unpack(flat_initial_params)

function build_gp(params)
    return GP(params.var_kernel * RationalKernel(α=params.α) ∘ ScaleTransform(params.λ))
end

# Define MLE objective
function objective(params)
    f = build_gp(unpack(params))
    return -logpdf(f(x, 0.001), y)
end

# Optimise using Optim
(best_f, best_θ, θ_initials, vistited_best_fs) = multi_start_hyperopt(objective, bfgs, 100, flat_initial_params)

# Extracting the optimal values of the parameters
optimal_params = unpack(best_θ)

# Final GP
f = build_gp(optimal_params)

# Finite projection of `f` at inputs `x`.
# Added Gaussian noise with variance 0.001.
fx = f(x, 0.001)

# Exact posterior given `y`. This is another GP.
p_fx = posterior(fx, y)

# Plot posterior.
plot!(plt_comp, collect(range_mul_factor) * 100 / max_load, p_fx(collect(range_mul_factor)); label="Posterior - " * string(RationalKernel))
```

### FBMKernel

```@example AbstractGPs

## Define GP prior with FBMKernel
flat_initial_params, unflatten = flatten((;
    var_kernel = bounded(10.0, 6.0, 1e8),
    λ = bounded(0.005, 0.0, 0.01),
    h = bounded(0.05, 0.0, 1.0),
))

# Construct a function to unpack flattened parameters and pull out the raw values.
unpack = ParameterHandling.value ∘ unflatten
params = unpack(flat_initial_params)

function build_gp(params)
    return GP(params.var_kernel * FBMKernel(h=params.h) ∘ ScaleTransform(params.λ))
end

# Define MLE objective
function objective(params)
    f = build_gp(unpack(params))
    return -logpdf(f(x, 0.001), y)
end

# Optimise using Optim
(best_f, best_θ, θ_initials, vistited_best_fs) = multi_start_hyperopt(objective, bfgs, 100, flat_initial_params)

# Extracting the optimal values of the parameters
optimal_params = unpack(best_θ)

# Final GP
f = build_gp(optimal_params)

# Finite projection of `f` at inputs `x`.
# Added Gaussian noise with variance 0.001.
fx = f(x, 0.001)

# Exact posterior given `y`. This is another GP.
p_fx = posterior(fx, y)

# Plot posterior.
plot!(plt_comp, collect(range_mul_factor) * 100 / max_load, p_fx(collect(range_mul_factor)); label="Posterior - " * string(FBMKernel))
```
