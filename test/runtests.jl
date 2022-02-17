using Clp
using OptimalBids
using OptimalBids.PowerModelsMarkets
using OptimalBids.PowerModelsMarkets: build_market, change_bids!, clear_market!, calculate_profit
using Test

import OptimalBids: Market, build_market, change_bids!, clear_market!, calculate_profit

DATA_DIR = joinpath(dirname(@__FILE__), "data")

include("core.jl")
include("powermodelsmarkets.jl")
