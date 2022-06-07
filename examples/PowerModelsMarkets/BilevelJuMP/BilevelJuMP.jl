using OptimalBids
using OptimalBids.PowerModelsMarkets
using JuMP

using BilevelJuMP
using Ipopt
using Gurobi
using QuadraticToBinary

using Plots # For some evaluation plots at the end
using Plots.PlotMeasures
using Downloads # To download Test Cases

# Needed fix for PowerModels interface
Base.getindex(v::BilevelJuMP.BilevelVariableRef, i::Int64) = v

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

# ### OptimalBids API

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
max_total_volume = 655.0
range_mul_factor = min_total_volume:1.0:max_total_volume
bid_range = [offer_weights .* [i] for i in range_mul_factor]
p_curve = profit_curve!(market, bid_range)

# Let's plot and see how the range profit evaluatiuon went:
plt_range = plot(collect(range_mul_factor) * 100 / max_load, p_curve,
    label="Regularized benchmark",
    ylabel="Profit (USD)",
    xlabel="Bid Volume (% Market Share)",
    legend=:outertopright,
    color="black",
    width=3,
    left_margin=10mm,
    bottom_margin=10mm,
    size=(900, 600)
);
plt_comp = deepcopy(plt_range);

# ### Bilevel NLP

# Make sure max bids are at their maximum
max_generations = offer_weights .* max_total_volume
change_bids!(market, max_generations)

# optimize
model = BilevelModel(Ipopt.Optimizer, 
    mode = BilevelJuMP.ProductMode(1e-5)
)

@variable(Upper(model), 0 <= qS[i=1:num_strategic_buses] <= max_generations[i], start = 2.0) # if `start=0.1` => false NLP infeastible

@constraint(Upper(model), sum(qS) <= max_total_volume)

pm = instantiate_model(market; jump_model=Lower(model))

@variable(Upper(model), -10_000 <= lambda[i=1:num_strategic_buses] <= 10_000, DualOf(sol(pm, 0, :bus, parse(Int64, bus_indexes[i]))[:lam_kcl_r]), start = 1)

gS = [PowerModels.var(pm, 0, :pg, parse(Int64, generator_indexes[i])) for i = 1:num_strategic_buses]

@constraint(Lower(model), gS .<= qS)

@objective(Upper(model), Max, - lambda'gS)

@elapsed optimize!(model)

# Final simulation
nlp_bids = value.(qS)
nlp_profit = profit_for_bid!(market, nlp_bids)

bid_range_nlp = [nlp_bids .* [i] for i in 0.1:0.1:10]
p_curve_nlp = profit_curve!(market, bid_range_nlp)
plot!(plt_comp, sum.(bid_range_nlp) * 100 / max_load, p_curve_nlp,
    label="Range Evaluation - NLP",
    color="purple"
);
scatter!(plt_comp, [sum(nlp_bids) * 100 / max_load], [nlp_profit],
    label="Bilevel Solution - NLP",
    color="purple"
)

# ### Bilevel Fortuny

# Make sure max bids are at their maximum
max_generations = offer_weights .* max_total_volume
change_bids!(market, max_generations)

# optimize
opt = QuadraticToBinary.Optimizer{Float64}(MOI.instantiate(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.2)))
model = BilevelModel(() -> opt, 
    mode = BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = 10000, dual_big_M = 10000)
)


@variable(Upper(model), 0 <= qS[i=1:num_strategic_buses] <= max_generations[i], start = nlp_bids[i])

@constraint(Upper(model), sum(qS) <= max_total_volume)

pm = instantiate_model(market; jump_model=Lower(model))

@variable(Upper(model), -10_000 <= lambda[i=1:num_strategic_buses] <= 10_000, DualOf(sol(pm, 0, :bus, parse(Int64, bus_indexes[i]))[:lam_kcl_r]), start = 1_000)

gS = [PowerModels.var(pm, 0, :pg, parse(Int64, generator_indexes[i])) for i = 1:num_strategic_buses]

@constraint(Lower(model), gS .<= qS)

@objective(Upper(model), Max, - lambda'gS)

@elapsed optimize!(model)

# Final simulation
fortuny_bids = value.(qS)
fortuny_profit = profit_for_bid!(market, fortuny_bids)

bid_range_fortuny = [fortuny_bids .* [i] for i in 0.1:0.1:10]
p_curve_fortuny = profit_curve!(market, bid_range_fortuny)
plot!(plt_comp, sum.(bid_range_fortuny) * 100 / max_load, p_curve_fortuny,
    label="Range Evaluation - Fortuny",
    color="green"
);
scatter!(plt_comp, [sum(fortuny_bids) * 100 / max_load], [fortuny_profit],
    label="Bilevel Solution - Fortuny",
    color="green"
)

plot!(plt_comp, xlim=(0,5))

# ### Bilevel SOS

# Make sure max bids are at their maximum
max_generations = offer_weights .* max_total_volume
change_bids!(market, max_generations)

# optimize
opt = QuadraticToBinary.Optimizer{Float64}(MOI.instantiate(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.2)))
model = BilevelModel(() -> opt, 
    mode = BilevelJuMP.SOS1Mode()
)

@variable(Upper(model), 0 <= qS[i=1:num_strategic_buses] <= max_generations[i], start = 0.05)

@constraint(Upper(model), sum(qS) <= max_total_volume)

pm = instantiate_model(market; jump_model=Lower(model))

@variable(Upper(model), -10_000 <= lambda[i=1:num_strategic_buses] <= 10_000, DualOf(sol(pm, 0, :bus, parse(Int64, bus_indexes[i]))[:lam_kcl_r]), start = 1_000)

gS = [PowerModels.var(pm, 0, :pg, parse(Int64, generator_indexes[i])) for i = 1:num_strategic_buses]

@constraint(Lower(model), gS .<= qS)

@objective(Upper(model), Max, - lambda'gS)

@elapsed optimize!(model)

# Final simulation
SOS_bids = value.(qS)
SOS_profit = profit_for_bid!(market, SOS_bids)

bid_range_SOS = [SOS_bids .* [i] for i in 0.1:0.1:10]
p_curve_SOS = profit_curve!(market, bid_range_SOS)
plot!(plt_comp, sum.(bid_range_SOS), p_curve_SOS,
    label="Range Evaluation - SOS",
    color="red"
);
scatter!(plt_comp, [sum(SOS_bids)], [SOS_profit],
    label="Bilevel Solution - SOS",
    color="red"
)