using Clp
using OptimalBids
using OptimalBids.PowerModelsMarkets
using Test
using Downloads
using ChainRulesCore
using ChainRulesTestUtils
using Nonconvex
using NonconvexIpopt # Nonconvex.@load Ipopt
using NonconvexBayesian # Nonconvex.@load BayesOpt

import OptimalBids: Market, build_market, change_bids!, clear_market!, calculate_profit

include("core.jl")
include("powermodelsmarkets.jl")
