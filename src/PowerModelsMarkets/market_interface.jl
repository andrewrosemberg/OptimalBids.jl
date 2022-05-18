
"""
PowerModelsMarket <: OptimalBids.Market

Energy-Market type that uses PowerModels' OPF to clear the auction.

Arguments:
- `network_data::Dict`: [PowerModels data structure](https://lanl-ansi.github.io/PowerModels.jl/stable/network-data/).
- `strategic_generators::Vector{NamedTuple{(:gen_index, :bus_index),Tuple{String,String}}}`: Vector of strategic generators' indexes and their bus indexes. 
- `result::Union{Dict,Missing}`: Market clearing result.
- `market_formulation`: Network formulation used in the PowerModels' auction clearing process (i.e. OPF).
- `opf_builder`: PowerModels opf builder.
- `solver`: JuMP optimization solver that should be able to solve the OPF created based on the passed `market_formulation`.
"""
mutable struct PowerModelsMarket <: OptimalBids.Market
    network_data::Dict
    strategic_generators::Vector{NamedTuple{(:gen_index, :bus_index),Tuple{String,String}}}
    result::Union{Dict,Missing}
    market_formulation
    opf_builder
    solver
end

"""
OptimalBids.build_market(
    ::Type{PowerModelsMarket},
    network_data,
    strategic_generators,
    solver;
    market_formulation=DCPPowerModel,
    opf_builder=PowerModels.build_opf,
    assert_consistency=true,
)

Creates Energy-Market of type `PowerModelsMarket`.

Arguments:
- `network_data::Dict`: [PowerModels data structure](https://lanl-ansi.github.io/PowerModels.jl/stable/network-data/).
- `strategic_generators::Vector{NamedTuple{(:gen_index, :bus_index),Tuple{String,String}}}`: Vector of strategic generators' indexes and their bus indexes. 
- `result::Union{Dict,Missing}`: Market clearing result.
- `market_formulation`: Network formulation used in the PowerModels' auction clearing process (i.e. OPF).
- `opf_builder`: PowerModels opf builder.
- `solver`: JuMP optimization solver that should be able to solve the OPF created based on the passed `market_formulation`.
- `assert_consistency=true`: Boolean to force check if strategic generators indexes (in `strategic_generators`) reference generators located at the specified `bus_index`.
"""
function OptimalBids.build_market(
    ::Type{PowerModelsMarket},
    network_data,
    strategic_generators,
    solver;
    market_formulation=DCPPowerModel,
    opf_builder=PowerModels.build_opf,
    assert_consistency=true,
)
    if assert_consistency && any(
        string(network_data["gen"][i.gen_index]["gen_bus"]) != i.bus_index for
        i in strategic_generators
    )
        throw(DomainError("gen_indexes & bus_indexes inconsistent"))
    end
    return PowerModelsMarket(
        network_data, strategic_generators, missing, market_formulation, opf_builder, solver
    )
end

"""
OptimalBids.build_market(
    market::Type{PowerModelsMarket},
    network_data::Dict,
    gen_indexes::AbstractVector{String},
    bus_indexes::AbstractVector{String},
    solver;
    market_formulation=DCPPowerModel,
    opf_builder=PowerModels.build_opf,
    kwards...,
)

Creates Energy-Market of type `PowerModelsMarket`.

Arguments:
- `network_data::Dict`: [PowerModels data structure](https://lanl-ansi.github.io/PowerModels.jl/stable/network-data/).
- `gen_indexes::AbstractVector{String}`: Vector of strategic generators' indexes.
- `bus_indexes::AbstractVector{String}`: Vector of strategic generators' bus indexes. 
- `result::Union{Dict,Missing}`: Market clearing result.
- `market_formulation`: Network formulation used in the PowerModels' auction clearing process (i.e. OPF).
- `opf_builder`: PowerModels opf builder.
- `solver`: JuMP optimization solver that should be able to solve the OPF created based on the passed `market_formulation`.
"""
function OptimalBids.build_market(
    market::Type{PowerModelsMarket},
    network_data::Dict,
    gen_indexes::AbstractVector{String},
    bus_indexes::AbstractVector{String},
    solver;
    market_formulation=DCPPowerModel,
    opf_builder=PowerModels.build_opf,
    kwards...,
)
    return build_market(
        market,
        network_data,
        [
            (; gen_index=gen_indexes[i], bus_index=bus_indexes[i]) for
            i in 1:length(gen_indexes)
        ],
        solver;
        market_formulation=market_formulation,
        opf_builder=opf_builder,
        kwards...,
    )
end

"""
OptimalBids.build_market(
    market::Type{PowerModelsMarket},
    network_data::Dict,
    gen_indexes::AbstractVector{String},
    solver;
    market_formulation=DCPPowerModel,
    opf_builder=PowerModels.build_opf,
)

Creates Energy-Market of type `PowerModelsMarket`.

Arguments:
- `network_data::Dict`: [PowerModels data structure](https://lanl-ansi.github.io/PowerModels.jl/stable/network-data/).
- `gen_indexes::AbstractVector{String}`: Vector of strategic generators' indexes.
- `result::Union{Dict,Missing}`: Market clearing result.
- `market_formulation`: Network formulation used in the PowerModels' auction clearing process (i.e. OPF).
- `opf_builder`: PowerModels opf builder.
- `solver`: JuMP optimization solver that should be able to solve the OPF created based on the passed `market_formulation`.
"""
function OptimalBids.build_market(
    market::Type{PowerModelsMarket},
    network_data::Dict,
    gen_indexes::AbstractVector{String},
    solver;
    market_formulation=DCPPowerModel,
    opf_builder=PowerModels.build_opf,
)
    return build_market(
        market,
        network_data,
        gen_indexes,
        [string(network_data["gen"][i]["gen_bus"]) for i in gen_indexes],
        solver;
        market_formulation=market_formulation,
        opf_builder=opf_builder,
        assert_consistency=false,
    )
end

function OptimalBids.change_bids!(market::PowerModelsMarket, new_bids::Vector{Float64})
    length(market.strategic_generators) != length(new_bids) &&
        throw(BoundsError("Number of generators & new_bids dont match"))
    for (i, generator) in enumerate(market.strategic_generators)
        market.network_data["gen"][generator.gen_index]["pmax"] = new_bids[i]
    end
    return nothing
end

"""
PowerModels.instantiate_model(market::PowerModelsMarket)

Instantiates PowerModels Model.
"""
function PowerModels.instantiate_model(market::PowerModelsMarket; jump_model=JuMP.Model())
    return instantiate_model(
        market.network_data,
        market.market_formulation,
        market.opf_builder;
        setting=Dict("output" => Dict("duals" => true)),
        jump_model=jump_model
    )
end

"""
OptimalBids.clear_market!(market::PowerModelsMarket)

Clears market and stores result's dictionary in `result`.
"""
function OptimalBids.clear_market!(market::PowerModelsMarket)
    market.result = optimize_model!(instantiate_model(market); optimizer=market.solver)
    return nothing
end

function OptimalBids.calculate_profit(market::PowerModelsMarket)
    cleared_volumes = map(market.strategic_generators) do generator
        return market.result["solution"]["gen"][generator.gen_index]["pg"]
    end
    clearing_prices = map(market.strategic_generators) do generator
        return -market.result["solution"]["bus"][generator.bus_index]["lam_kcl_r"]
    end
    return (;
        cleared_volumes=cleared_volumes,
        clearing_prices=clearing_prices,
        profit=cleared_volumes .* clearing_prices,
    )
end
