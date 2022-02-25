module OptimalBids

export profit_for_bid!, profit_curve!

"""
    Market

Mutable Structure of a market instance.
"""
abstract type Market end

"""
    build_market(::Type{Market}, params...) -> Market

Builds market of type Market using provided parameters (`params`).
"""
function build_market(typ::Type{M}, params...) where {M<:Market}
    throw(MethodError(build_market, (typ, params...)))
end

"""
    change_bids!(market::Market, new_bids::Vector)

Changes strategic agent's bids in the market to `new_bids`.
"""
function change_bids!(market::M, new_bids::Vector) where {M<:Market}
    throw(MethodError(change_bids!, (market, new_bids)))
end

"""
    clear_market!(market::Market)

Clears the market.
"""
function clear_market!(market::M) where {M<:Market}
    throw(MethodError(clear_market!, (market)))
end

"""
    calculate_profit(market::Market) -> NamedTuple{(:cleared_volumes, :clearing_prices, :profit), Tuple{Vector{Int}, Vector{Int}, Vector{Int}}}

Retrieves strategic agent's cleared volumes and prices from the market and calculates per bid profit.
"""
function calculate_profit(market::M) where {M<:Market}
    throw(MethodError(calculate_profit, (market)))
end

"""
    profit_for_bid!(market::Market, new_bids::Vector{Any}) -> Float64

Calculates overall profit when market is cleared with `new_bids`. This function will sequentiall call `change_bids!`, `clear_market!` and calculate_profit.
"""
function profit_for_bid!(market::Market, new_bids::Vector{<:Any})
    change_bids!(market, new_bids)
    clear_market!(market)
    return sum(calculate_profit(market).profit)
end

"""
    profit_curve!(market::Market, range_new_bids::Vector{Vector{Any}}) -> Vector{Float64}

Constructs profit curve for bids provided in `range_new_bids`.
"""
function profit_curve!(market::Market, range_new_bids::Vector{<:AbstractVector})
    map(range_new_bids) do bids
        return profit_for_bid!(market, bids)
    end
end

# Include SubModules
include("PowerModelsMarkets/PowerModelsMarkets.jl")

end
