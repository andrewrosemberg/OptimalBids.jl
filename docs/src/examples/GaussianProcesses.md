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
max_total_volume = 555.0
range_mul_factor = min_total_volume:1.0:max_total_volume
bid_range = [offer_weights .* [i] for i in range_mul_factor]
p_curve = profit_curve!(market, bid_range)

```

## AbstractGPs Interface

### Data

```@example AbstractGPs

# Randomly choose observed data
idx = rand(1:size(range_mul_factor,1), 10)
x = collect(range_mul_factor)[idx]
y = p_curve[idx]

plt_comp = plot(collect(range_mul_factor), p_curve,
    label="Range Evaluation",
    ylabel="Profit (\$)",
    xlabel="Multiplicative Factor",
    legend=:outertopright,
    color="black",
    width=3,
);
scatter!(plt_comp, x, y; label="Data")

```

### Matern32Kernel

```@example AbstractGPs

## Define GP prior with Matern32Kernel
flat_initial_params, unflatten = flatten((;
    var_kernel = positive(0.6),
    λ = positive(50000.5),
))

# Construct a function to unpack flattened parameters and pull out the raw values.
unpack = ParameterHandling.value ∘ unflatten
params = unpack(flat_initial_params)

function build_gp(params)
    return GP(params.var_kernel * Matern52Kernel() ∘ ScaleTransform(params.λ))
end

# Define MLE objective
function objective(params)
    f = build_gp(params)
    return -logpdf(f(x, 0.001), y)
end

# Optimise using Optim
training_results = Optim.optimize(
    objective ∘ unpack,
    θ -> only(Zygote.gradient(objective ∘ unpack, θ)),
    [log(0.001); log(1.0)],
    [log(3.0); log(100000000.0)],
    flat_initial_params,
    Fminbox(BFGS(
        alphaguess = Optim.LineSearches.InitialStatic(scaled=true),
        linesearch = Optim.LineSearches.BackTracking(),
    )),
    Optim.Options(show_trace = true);
    inplace=false,
)

# Extracting the optimal values of the parameters
optimal_params = unpack(training_results.minimizer)

# Final GP
f = build_gp(optimal_params)

# Finite projection of `f` at inputs `x`.
# Added Gaussian noise with variance 0.001.
fx = f(x, 0.001)

# Exact posterior given `y`. This is another GP.
p_fx = posterior(fx, y)

# Plot posterior.
plot!(plt_comp, range_mul_factor, p_fx; label="Posterior - " * string(Matern32Kernel))

```