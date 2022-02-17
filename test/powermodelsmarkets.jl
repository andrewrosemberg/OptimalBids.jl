@testset "PowerModelsMarkets" begin
    @testset "Case: $case_name" for case_name in ["case118.m", "case5_re_uc.m"]
        case_file_path = joinpath(DATA_DIR, case_name)

        # Choice of number of strategic bidding locations
        percentage_buses = 0.1

        # Read data
        network_data = PowerModels.parse_file(case_file_path)

        # get indexes
        bus_indexes = collect(keys(network_data["bus"]))
        num_buses = length(bus_indexes)
        num_strategic_buses = ceil(Int, percentage_buses * num_buses)
        first_index_bus = floor(Int, (num_buses - num_strategic_buses) / 2)
        bus_indexes = bus_indexes[first_index_bus:(first_index_bus + num_strategic_buses - 1)] # lets grab the middle generators

        # add new strategic generators
        generator_indexes = [
            add_generator(network_data, parse(Int, bus_idx)) for bus_idx in bus_indexes
        ]

        # test PowerModelsMarket functions
        market = nothing
        @testset "build_market" begin
            @test_throws MethodError build_market(
                PowerModelsMarket,
                network_data,
                collect(1:num_strategic_buses),
                collect(1:num_strategic_buses),
                Clp.Optimizer,
            )
            @test_throws BoundsError build_market(
                PowerModelsMarket,
                network_data,
                generator_indexes,
                bus_indexes[1:(end - 1)],
                Clp.Optimizer,
            )
            @test_throws DomainError build_market(
                PowerModelsMarket,
                network_data,
                generator_indexes,
                [bus_indexes[1:(end - 1)]; "44"],
                Clp.Optimizer,
            )

            sg_aux =
                build_market(
                    PowerModelsMarket,
                    network_data,
                    generator_indexes,
                    bus_indexes,
                    Clp.Optimizer,
                ).strategic_generators

            market = build_market(
                PowerModelsMarket,
                network_data,
                generator_indexes,
                bus_indexes,
                Clp.Optimizer,
            )
            @test all(market.strategic_generators .== sg_aux)
        end

        initial_bids = collect(0.01:0.01:(num_strategic_buses * 0.01))
        @testset "change_bids!" begin
            @test_throws BoundsError change_bids!(
                market, collect(1.0:(num_strategic_buses - 1))
            )

            change_bids!(market, initial_bids)
            @test [
                market.network_data["gen"][i.gen_index]["pmax"] for
                i in market.strategic_generators
            ] == initial_bids
        end

        @testset "clear_market!" begin
            @test ismissing(market.result)
            clear_market!(market)
            @test !ismissing(market.result)
        end

        @testset "calculate_profit" begin
            @test length(calculate_profit(market).profit) == num_strategic_buses
        end

        @testset "profit_for_bid!" begin
            @test profit_for_bid!(market, zeros(num_strategic_buses)) == 0.0
        end

        @testset "profit_curve!" begin
            bid_range = [initial_bids .* [i] for i in 0.0:10]
            p_curve = profit_curve!(market, bid_range)
            @test length(p_curve) == length(bid_range)
            @test maximum(p_curve) > minimum(p_curve)
        end
    end
end
