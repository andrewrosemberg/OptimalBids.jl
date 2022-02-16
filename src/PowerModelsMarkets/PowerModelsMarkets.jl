"""
    PowerModelsMarkets

SubModule that implements interfact with PowerModels.jl
"""
module PowerModelsMarkets
using PowerModels
using JuMP

import ..OptimalBids: Market, build_market, change_bids!, clear_market!, calculate_profit

export PowerModelsMarket,
    build_market, change_bids!, clear_market!, calculate_profit, add_generator

include("case_modifiers.jl")

mutable struct PowerModelsMarket <: OptimalBids.Market
    network_data::Dict
    strategic_generators::Vector{NamedTuple{(:gen_index, :bus_index),Tuple{Int64,Int64}}}
    result::Union{Dict,Missing}
    market_formulation::AbstractActivePowerModel
    opf_builder::Fuction
    solver
end

function OptimalBids.build_market(
    ::Type{PowerModelsMarket},
    network_data,
    strategic_generators,
    solver;
    market_formulation=DCPPowerModel,
    opf_builder=PowerModels.build_opf,
)
    return PowerModelsMarket(
        network_data, strategic_generators, missing, market_formulation, opf_builder, solver
    )
end

function OptimalBids.change_bids!(market::PowerModelsMarket, bids::Vector{Float64})
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
        return result["solution"]["gen"][generator.gen_index]["pg"]
    end
    clearing_price = map(market.strategic_generators) do generator
        return -result["solution"]["bus"][generator.bus_index]["lam_kcl_r"]
    end
    return (;
        cleared_volumes=cleared_volumes,
        clearing_price=clearing_price,
        profit=market.cleared_volumes .* market.cleared_prices,
    )
end
end
