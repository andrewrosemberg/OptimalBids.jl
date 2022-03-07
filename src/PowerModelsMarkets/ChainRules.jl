function ChainRulesCore.rrule(::typeof(profit_for_bid!), market::PowerModelsMarket, inputs::Vector{<:Any})
    change_bids!(market, inputs)
    clear_market!(market)
    clearing = calculate_profit(market)

    y = sum(clearing.profit)
    
    function profit_for_bid_pullback(ȳ)
        īnputs = deepcopy(clearing.clearing_prices)
        for i in 1:length(inputs)
            if isapprox(clearing.cleared_volumes[i], inputs[i], rtol=0.1)
                īnputs[i] *= ȳ
            else
                īnputs[i] = 0.0
            end
        end
        m̄arket = @not_implemented("derivative w.r.t market")
        return NoTangent(), m̄arket, īnputs
    end
    
    return y, profit_for_bid_pullback
end