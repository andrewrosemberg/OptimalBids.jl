"""
    analytical_gradient(::typeof(profit_for_bid!), clearing::NamedTuple, inputs::Vector{<:Any}) -> Vector{Float64}

Analytical gradient of `profit_for_bid!` w.r.t inputs, knowing clearing.
"""
function analytical_gradient(::typeof(profit_for_bid!), clearing::NamedTuple, inputs::Vector{<:Any})
    Δ = deepcopy(clearing.clearing_prices)
    for i in 1:length(inputs)
        if !isapprox(clearing.cleared_volumes[i], inputs[i], rtol=0.1)
            Δ[i] = 0.0
        end
    end
    return Δ
end

"""
    profit_and_gradient_for_bid!(market::Market, new_bids::Vector{Any}) -> Float64, Float64

Calculates overall profit when market is cleared with `new_bids` and its gradient w.r.t the vector of bids.
"""
function profit_and_gradient_for_bid!(market::Market, new_bids::Vector{<:Any})
    change_bids!(market, new_bids)
    clear_market!(market)
    clearing = calculate_profit(market)

    y = sum(clearing.profit)
    Δ = analytical_gradient(profit_for_bid!, clearing, new_bids)

    return y, Δ
end

function ChainRulesCore.rrule(f::typeof(profit_for_bid!), market::PowerModelsMarket, inputs::Vector{<:Any})
    y, Δ = profit_and_gradient_for_bid!(market, inputs)
    
    function profit_for_bid_pullback(ȳ)
        īnputs = Δ * ȳ
        m̄arket = @not_implemented("derivative w.r.t market")
        return NoTangent(), m̄arket, īnputs
    end
    
    return y, profit_for_bid_pullback
end