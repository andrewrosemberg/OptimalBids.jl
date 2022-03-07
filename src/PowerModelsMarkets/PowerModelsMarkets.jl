"""
    PowerModelsMarkets

SubModule that implements interfact with PowerModels.jl
"""
module PowerModelsMarkets
using ..OptimalBids
using PowerModels
using JuMP
using ChainRulesCore

import ..OptimalBids: Market, build_market, change_bids!, clear_market!, calculate_profit
using Reexport: Reexport

export PowerModelsMarket,
    build_market, change_bids!, clear_market!, calculate_profit, add_generator

Reexport.@reexport using PowerModels

include("case_modifiers.jl")
include("market_interface.jl")
include("ChainRules.jl")

end
