module OptimalBids

abstract type Market end

"""
    build_market(::Type{Market}, params...) -> Market

Builds market of type Market using provided parameters (`params`).
"""
function build_market(::Type{Market}, params...)
    @error "Market builder not defined for a market of type $Market."
end

"""
    change_bids!(market::Market, new_bids::Vector{Any})

Changes strategic agent's bids in the market to `new_bids`.
"""
function change_bids!(market::Market, new_bids::Vector{Any})
    @error "function change_bids! not defined for a market of type $Market."
end

"""
    clear_market!(market::Market)

Clears the market.
"""
function clear_market!(market::Market)
    @error "function clear_market! not defined for a market of type $Market."
end

"""
    calculate_profit(market::Market) -> NamedTuple{(:cleared_volumes, :clearing_price, :profit), Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}}

Retrieves strategic agent's cleared volumes and prices from the market and calculates per bid profit.
"""
function calculate_profit(market::Market)
    @error "function calculate_profit not defined for market of type $Market."
end

"""
    profit_for_bid!(market::Market, new_bids::Vector{Any}) -> Float64

Calculates overall profit when market is cleared with `new_bids`. This function will sequentiall call `change_bids!`, `clear_market!` and calculate_profit.
"""
function profit_for_bid!(market::Market, new_bids::Vector{Any})
    change_bids!(market, new_bids)
    clear_market!(market)
    return sum(calculate_profit(market).profit)
end

"""
    profit_curve!(market::Market, range_new_bids::Vector{Vector{Any}}) -> Vector{Float64}

Constructs profit curve for bids provided in `range_new_bids`.
"""
function profit_curve!(market::Market, range_new_bids::Vector{Vector{Any}})
    map(range_new_bids) do bids
        return profit_for_bid!(market, bids)
    end
end
