using OptimalBids
using Test

import OptimalBids: Market, build_market, change_bids!, clear_market!, calculate_profit, profit_for_bid!, profit_curve!

include("core.jl")
