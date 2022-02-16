@testset "Core OptimalBids.jl" begin
    mutable struct MyMarket <: OptimalBids.Market
        bids::Vector{Float64}
        cleared_volumes::Union{Vector{Float64}, Missing}
        cleared_prices::Union{Vector{Float64}, Missing}
    end

    initial_bids = collect(1.0:10.0)
    market = MyMarket(initial_bids, missing, missing)
    @testset "build_market" begin
        @test_throws MethodError OptimalBids.build_market(MyMarket, initial_bids)

        function OptimalBids.build_market(::Type{MyMarket}, bids::Vector{Float64})
            MyMarket(bids, missing, missing)
        end

        market = OptimalBids.build_market(MyMarket, initial_bids) 
    end

    @testset "change_bids!" begin
        new_bids = collect(5.0:15.0)
        @test_throws MethodError OptimalBids.change_bids!(MyMarket, new_bids)

        function OptimalBids.change_bids!(market::MyMarket, bids::Vector{Float64})
            market.bids = bids
        end

        OptimalBids.change_bids!(market, new_bids) 
    end
    
    @testset "clear_market!" begin
        @test_throws MethodError OptimalBids.clear_market!(MyMarket)

        function OptimalBids.clear_market!(market::MyMarket)
            market.cleared_volumes = market.bids ./ 2
            market.cleared_prices = ones(length(market.bids))
        end

        OptimalBids.clear_market!(market) 
    end
    
    @testset "calculate_profit" begin
        @test_throws MethodError OptimalBids.calculate_profit(MyMarket)

        function OptimalBids.calculate_profit(market::MyMarket)
            (; cleared_volumes = market.cleared_volumes, clearing_price = market.cleared_prices, profit = sum(market.cleared_volumes .* market.cleared_prices))
        end

        OptimalBids.calculate_profit(market) 
    end
    
    @testset "profit_for_bid!" begin
        @test OptimalBids.profit_for_bid!(market, initial_bids) == sum(initial_bids) / 2
    end

    @testset "profit_curve!" begin
        @test OptimalBids.profit_curve!(market, [initial_bids .* [i] for i in 1:10]) == [sum(initial_bids) * i / 2 for i in 1:10]
    end
end
