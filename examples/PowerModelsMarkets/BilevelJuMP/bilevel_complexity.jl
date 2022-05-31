using OptimalBids
using OptimalBids.PowerModelsMarkets
using JuMP
using JuMP.Containers

using BilevelJuMP
using Ipopt
using Gurobi
using QuadraticToBinary

using Plots # For some evaluation plots at the end
using Plots.PlotMeasures
using Downloads # To download Test Cases
using Random

# Needed fix for PowerModels interface
Base.getindex(v::BilevelJuMP.BilevelVariableRef, i::Int64) = v

# ### Case Definition
cases = Dict(
    5 => "pglib_opf_case5_pjm.m",
    # 14 => "pglib_opf_case14_ieee.m",
    30 => "pglib_opf_case30_ieee.m",
    # 73 => "pglib_opf_case73_ieee_rts.m",
    # 118 => "pglib_opf_case118_ieee.m",
)

time_solve = DenseAxisArray(zeros(length(keys(cases)), 3), collect(keys(cases)), [:NLP, :FORTUNY, :SOS1])
obj_val = DenseAxisArray(zeros(length(keys(cases)), 3), collect(keys(cases)), [:NLP, :FORTUNY, :SOS1])

# Read market data from IEEE 118 bus case
DATA_DIR = mktempdir()
env = Gurobi.Env()
rng = MersenneTwister(666)

for case_name in values(cases)
    case_file_path = joinpath(DATA_DIR, case_name)
    if !isfile(case_file_path)
        case_file_path = joinpath(DATA_DIR, case_name)
        Downloads.download("https://raw.githubusercontent.com/power-grid-lib/pglib-opf/01681386d084d8bd03b429abcd1ee6966f68b9a3/" * case_name, case_file_path)
    end
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
    bus_indexes = rand(rng, bus_indexes, num_strategic_buses)
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
        optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0),
    )

    min_total_volume = 0.0
    max_total_volume = 655.0
    range_mul_factor = min_total_volume:1.0:max_total_volume

    # ### Bilevel NLP

    # Make sure max bids are at their maximum
    max_generations = ones(num_strategic_buses) * max_total_volume
    change_bids!(market, max_generations)

    # optimize
    model = BilevelModel(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0), 
        mode = BilevelJuMP.ProductMode(1e-5)
    )

    @variable(Upper(model), 0 <= qS[i=1:num_strategic_buses] <= max_generations[i], start = 0.1) # if `start=0.1` => false NLP infeastible

    @constraint(Upper(model), sum(qS) <= max_total_volume)

    pm = instantiate_model(market; jump_model=Lower(model))

    @variable(Upper(model), -10_000 <= lambda[i=1:num_strategic_buses] <= 10_000, DualOf(sol(pm, 0, :bus, parse(Int64, bus_indexes[i]))[:lam_kcl_r]), start = 1)

    gS = [var(pm, 0, :pg, parse(Int64, generator_indexes[i])) for i = 1:num_strategic_buses]

    @constraint(Lower(model), gS .<= qS)

    @objective(Upper(model), Max, - lambda'gS)

    time_solve[num_buses, :NLP] = @elapsed optimize!(model)

    obj_val[num_buses, :NLP] = profit_for_bid!(market, value.(qS))

    finalize(backend(model).optimizer.model)

    # ### Bilevel Fortuny

    # Make sure max bids are at their maximum
    max_generations = ones(num_strategic_buses) * max_total_volume
    change_bids!(market, max_generations)

    # optimize
    opt = QuadraticToBinary.Optimizer{Float64}(MOI.instantiate(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "MIPGap" => 0.05)))
    model = BilevelModel(() -> opt, 
        mode = BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = 10000, dual_big_M = 10000)
    )


    @variable(Upper(model), 0 <= qS[i=1:num_strategic_buses] <= max_generations[i], start = 0.005)

    @constraint(Upper(model), sum(qS) <= max_total_volume)

    pm = instantiate_model(market; jump_model=Lower(model))

    @variable(Upper(model), -10_000 <= lambda[i=1:num_strategic_buses] <= 10_000, DualOf(sol(pm, 0, :bus, parse(Int64, bus_indexes[i]))[:lam_kcl_r]), start = 1_000)

    gS = [var(pm, 0, :pg, parse(Int64, generator_indexes[i])) for i = 1:num_strategic_buses]

    @constraint(Lower(model), gS .<= qS)

    @objective(Upper(model), Max, - lambda'gS)

    time_solve[num_buses, :FORTUNY] = @elapsed optimize!(model)
    obj_val[num_buses, :FORTUNY] = profit_for_bid!(market, value.(qS))

    # ### Bilevel SOS

    # Make sure max bids are at their maximum
    max_generations = ones(num_strategic_buses) * max_total_volume
    change_bids!(market, max_generations)

    # optimize
    opt = QuadraticToBinary.Optimizer{Float64}(MOI.instantiate(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "MIPGap" => 0.05)))
    model = BilevelModel(() -> opt, 
        mode = BilevelJuMP.SOS1Mode()
    )

    @variable(Upper(model), 0 <= qS[i=1:num_strategic_buses] <= max_generations[i], start = 0.05)

    @constraint(Upper(model), sum(qS) <= max_total_volume)

    pm = instantiate_model(market; jump_model=Lower(model))

    @variable(Upper(model), -10_000 <= lambda[i=1:num_strategic_buses] <= 10_000, DualOf(sol(pm, 0, :bus, parse(Int64, bus_indexes[i]))[:lam_kcl_r]), start = 1_000)

    gS = [var(pm, 0, :pg, parse(Int64, generator_indexes[i])) for i = 1:num_strategic_buses]

    @constraint(Lower(model), gS .<= qS)

    @objective(Upper(model), Max, - lambda'gS)

    time_solve[num_buses, :SOS1] = @elapsed optimize!(model)
    obj_val[num_buses, :SOS1] = profit_for_bid!(market, value.(qS))

    println("SOLVE TIME")
    println(time_solve)
    println("OBJ VALUE")
    println(obj_val)
end

finalize(env)
