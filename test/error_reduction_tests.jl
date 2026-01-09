using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase, Test

@testset "Error Reduction Function Tests" begin
    # Define a multi-component problem (Lotka-Volterra)
    function lotka(du, u, p, t)
        du[1] = 1.5 * u[1] - u[1] * u[2]
        du[2] = -3 * u[2] + u[1] * u[2]
        return du[3] = u[1] + u[2]  # Third component is sum of first two
    end
    prob_lotka = ODEProblem(lotka, [1.0, 1.0, 2.0], (0.0, 1.0))

    # Get reference solution
    ref_sol = solve(prob_lotka, Vern9(); abstol = 1.0e-14, reltol = 1.0e-14)
    test_sol = TestSolution(ref_sol)

    # Get a solution with some error
    sol = solve(prob_lotka, DP5(); abstol = 1.0e-4, reltol = 1.0e-4)

    @testset "appxtrue with default_reduction" begin
        # Test that default_reduction is identity
        @test DiffEqDevTools.default_reduction([1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]
        @test DiffEqDevTools.default_reduction(5.0) == 5.0

        # Test appxtrue with default reduction (should match original behavior)
        errsol_default = appxtrue(sol, ref_sol)
        @test haskey(errsol_default.errors, :final)
        @test haskey(errsol_default.errors, :l2)
        @test haskey(errsol_default.errors, :l∞)
    end

    @testset "appxtrue with TestSolution and reduction" begin
        # Test appxtrue with TestSolution
        errsol_default = appxtrue(sol, test_sol)

        # Test with reduction that selects only the third component
        reduction_third = x -> x[3]
        errsol_third = appxtrue(sol, test_sol; reduction = reduction_third)

        # Errors should be different (unless third component happens to have same error)
        @test haskey(errsol_third.errors, :final)
        @test haskey(errsol_third.errors, :l2)
        @test haskey(errsol_third.errors, :l∞)

        # Test with reduction that selects first component
        reduction_first = x -> x[1]
        errsol_first = appxtrue(sol, test_sol; reduction = reduction_first)
        @test haskey(errsol_first.errors, :final)
    end

    @testset "appxtrue with AbstractODESolution and reduction" begin
        # Test appxtrue with another ODESolution
        errsol_default = appxtrue(sol, ref_sol)

        # Test with reduction that computes norm of first two components
        reduction_first_two = x -> x[1:2]
        errsol_first_two = appxtrue(sol, ref_sol; reduction = reduction_first_two)
        @test haskey(errsol_first_two.errors, :final)
        @test haskey(errsol_first_two.errors, :l2)
    end

    @testset "WorkPrecision with reduction" begin
        abstols = 1 ./ 10 .^ (4:6)
        reltols = 1 ./ 10 .^ (1:3)

        # Test WorkPrecision with default reduction
        wp_default = WorkPrecision(
            prob_lotka, DP5(), abstols, reltols;
            appxsol = test_sol, name = "DP5-default"
        )
        @test wp_default.name == "DP5-default"

        # Test WorkPrecision with custom reduction (third component only)
        reduction_third = x -> x[3]
        wp_third = WorkPrecision(
            prob_lotka, DP5(), abstols, reltols;
            appxsol = test_sol, name = "DP5-third", reduction = reduction_third
        )
        @test wp_third.name == "DP5-third"

        # Verify both produced valid errors
        @test all(!isnan, wp_default.errors.final)
        @test all(!isnan, wp_third.errors.final)
    end

    @testset "WorkPrecisionSet with reduction" begin
        abstols = 1 ./ 10 .^ (4:6)
        reltols = 1 ./ 10 .^ (1:3)

        setups = [Dict(:alg => DP5()), Dict(:alg => Tsit5())]

        # Test WorkPrecisionSet with default reduction
        wp_set_default = WorkPrecisionSet(
            prob_lotka, abstols, reltols, setups;
            appxsol = test_sol
        )
        @test wp_set_default.N == 2

        # Test WorkPrecisionSet with custom reduction
        reduction_first = x -> x[1]
        wp_set_first = WorkPrecisionSet(
            prob_lotka, abstols, reltols, setups;
            appxsol = test_sol, reduction = reduction_first
        )
        @test wp_set_first.N == 2
    end

    @testset "Shootout with reduction" begin
        setups = [Dict(:alg => DP5()), Dict(:alg => Tsit5())]

        # Test Shootout with default reduction
        shoot_default = Shootout(prob_lotka, setups; appxsol = test_sol)
        @test shoot_default.N == 2

        # Test Shootout with custom reduction (second component only)
        reduction_second = x -> x[2]
        shoot_second = Shootout(
            prob_lotka, setups;
            appxsol = test_sol, reduction = reduction_second
        )
        @test shoot_second.N == 2

        # Both should have valid errors
        @test all(!isnan, shoot_default.errors)
        @test all(!isnan, shoot_second.errors)
    end

    @testset "Reduction function validation" begin
        # Test various reduction functions work correctly
        diff_vec = [0.1, 0.2, 0.3]

        # Identity reduction
        @test DiffEqDevTools.default_reduction(diff_vec) == diff_vec

        # Single element selection
        @test (x -> x[1])(diff_vec) == 0.1
        @test (x -> x[2])(diff_vec) == 0.2
        @test (x -> x[3])(diff_vec) == 0.3

        # Slice selection
        @test (x -> x[1:2])(diff_vec) == [0.1, 0.2]
        @test (x -> x[2:3])(diff_vec) == [0.2, 0.3]

        # Custom function of the solution
        @test (x -> sum(x))(diff_vec) ≈ 0.6
        @test (x -> x[1] - x[2])(diff_vec) ≈ -0.1
    end
end
