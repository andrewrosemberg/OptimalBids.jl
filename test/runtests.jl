using Clp
using OptimalBids
using OptimalBids.PowerModelsMarkets
using Test
using Nonconvex
using NonconvexIpopt # Nonconvex.@load Ipopt
using NonconvexBayesian # Nonconvex.@load BayesOpt

import OptimalBids: Market, build_market, change_bids!, clear_market!, calculate_profit

DATA_DIR = joinpath(dirname(@__FILE__), "data")

include("core.jl")
include("powermodelsmarkets.jl")
