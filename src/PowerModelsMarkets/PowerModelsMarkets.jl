"""
    PowerModelsMarkets

SubModule that implements interfact with PowerModels.jl
"""
module PowerModelsMarkets
using ..OptimalBids
using PowerModels
using JuMP

import ..OptimalBids: Market, build_market, change_bids!, clear_market!, calculate_profit
using Reexport: Reexport

export PowerModelsMarket,
    build_market, change_bids!, clear_market!, calculate_profit, add_generator

Reexport.@reexport using PowerModels

include("case_modifiers.jl")

mutable struct PowerModelsMarket <: OptimalBids.Market
    network_data::Dict
    strategic_generators::Vector{NamedTuple{(:gen_index, :bus_index),Tuple{String,String}}}
    result::Union{Dict,Missing}
    market_formulation
    opf_builder
    solver
end

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

function OptimalBids.change_bids!(market::PowerModelsMarket, bids::Vector{Float64})
    length(market.strategic_generators) != length(bids) &&
        throw(BoundsError("Number of generators & bids dont match"))
    for (i, generator) in enumerate(market.strategic_generators)
        market.network_data["gen"][generator.gen_index]["pmax"] = bids[i]
    end
    return nothing
end

function OptimalBids.clear_market!(market::PowerModelsMarket)
    pm = PowerModels.instantiate_model(
        market.network_data,
        market.market_formulation,
        market.opf_builder;
        setting=Dict("output" => Dict("duals" => true)),
    )
    market.result = optimize_model!(pm; optimizer=market.solver)
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
end
