@testset "PowerModelsMarkets" begin
    @testset "Case: $case_name" for case_name in ["pglib_opf_case118_ieee.m", "pglib_opf_case5_pjm.m"]
        # Download test case
        DATA_DIR = mktempdir()
        case_file_path = joinpath(DATA_DIR, case_name)
        Downloads.download(joinpath("https://raw.githubusercontent.com/power-grid-lib/pglib-opf/01681386d084d8bd03b429abcd1ee6966f68b9a3", case_name), case_file_path)

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

        @testset "rrule profit_for_bid!" begin
            test_rrule(profit_for_bid!, market ⊢ ChainRulesCore.NoTangent(), initial_bids)
            test_rrule(profit_for_bid!, market ⊢ ChainRulesCore.NoTangent(), initial_bids * 1e7)
        end

        min_total_volume = 0.0
        max_total_volume = 100.0
        best_profit = 0.0
        @testset "profit_curve!" begin
            bid_range = [initial_bids .* [i] for i in min_total_volume:max_total_volume]
            p_curve = profit_curve!(market, bid_range)
            best_profit = maximum(p_curve)
            @test length(p_curve) == length(bid_range)
            @test best_profit > minimum(p_curve)
        end

        @testset "Profit Maximization using Nonconvex" begin
            function profit_function(total_volume)
                initial_bids = collect(0.01:0.01:(num_strategic_buses * 0.01))
                return - profit_for_bid!(market, initial_bids .* total_volume[1])
            end

            # Max Number of Iterations for the solution method
            maxiter = 10

            model = Model()
            set_objective!(model, profit_function, flags = [:expensive])
            addvar!(model, [min_total_volume], [max_total_volume])
            add_ineq_constraint!(model, x -> -1) # Errors when no inequality is added! 

            alg = BayesOptAlg(IpoptAlg())
            options = BayesOptOptions(
                sub_options = IpoptOptions(),
                maxiter = maxiter, ftol = 1e-4, ctol = 1e-5,
            )
            r = optimize(model, alg, [min_total_volume], options = options)

            @test r.minimizer[1] >= min_total_volume && r.minimizer[1] <= max_total_volume
            # Assumption that the range evaluation found something at most 1.5x lower than the true maximum
            @test -r.minimum >= 0.0 && -r.minimum <= best_profit * 1.5
            @test r.niters <= maxiter
        end
    end
end
