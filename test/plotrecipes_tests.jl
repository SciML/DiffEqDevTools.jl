using Test
using OrdinaryDiffEq, StochasticDiffEq, DiffEqDevTools, Plots
import SDEProblemLibrary: prob_sde_additivesystem

using Random
Random.seed!(123)

gr()

@testset "ODE WorkPrecisionSet" begin
    # Linear ODE system
    f_rode_lin = function (du, u, p, t)
        @inbounds for i in eachindex(u)
            du[i] = cos(t) * u[i]
        end
    end

    f_rode_lin_analytic = function (u₀, p, t)
        u₀ * exp(sin(t))
    end

    tspan = (0.0, 10.0)
    prob = ODEProblem(
        ODEFunction(f_rode_lin, analytic = f_rode_lin_analytic), rand(10, 10),
        tspan)

    abstols = 1.0 ./ 10.0 .^ (3:8)
    reltols = 1.0 ./ 10.0 .^ (0:5)

    setups = [Dict(:alg => DP5())
              Dict(:alg => Tsit5())]
    wp_names = ["DP5", "Tsit5"]
    wp = WorkPrecisionSet(prob, abstols, reltols, setups, names = wp_names,
        save_everystep = false,
        numruns = 100)

    plt = @test_nowarn plot(wp)
    @test plt[1][1][:x] ≈ getproperty(wp[1].errors, wp[1].error_estimate)
    @test plt[1][2][:x] ≈ getproperty(wp[2].errors, wp[2].error_estimate)
    @test plt[1][1][:label] == wp_names[1]
    @test plt[1][2][:label] == wp_names[2]
    @test_nowarn plot(wp, color = [1 2])
    @test_nowarn plot(wp, color = :blue)
    @test_throws ArgumentError plot(wp, view = :dt_convergence)

    @test_nowarn plot(wp, x = :abstols, y = :final)
    @test_nowarn plot(wp, x = :reltols, y = :nf)
    @test_nowarn plot(wp, x = :naccept, y = :nreject)
    @test_throws ArgumentError plot(wp, x = :notakey, y = :final)

    dts = 1.0 ./ 2.0 .^ ((1:length(reltols)) .+ 1)
    setups = [Dict(:alg => Euler(), :dts => dts)
              Dict(:alg => Heun(), :dts => dts)
              Dict(:alg => Tsit5(), :dts => dts, :adaptive => false)
              Dict(:alg => Tsit5())]
    wp_names = ["Euler", "Heun", "Tsit5 fixed step", "Tsit5 adaptive"]
    wp = WorkPrecisionSet(prob, abstols, reltols, setups, names = wp_names,
        save_everystep = false,
        numruns = 100)

    plt = @test_nowarn plot(wp)
    @test all(plt[1][i][:x] ≈ getproperty(wp[i].errors, wp[i].error_estimate) for i in 1:4)
    @test all(plt[1][i][:label] == wp_names[i] for i in 1:4)

    @test_nowarn plot(wp, x = :dts, y = :final)

    plt = @test_nowarn plot(wp, view = :dt_convergence, legend = :bottomright)
    @test all(plt[1][i][:x] == plt[1][i + 3][:x] == dts == wp.setups[i][:dts] for i in 1:3)
    @test all(plt[1][i + 3][:y] ≈ getproperty(wp[i].errors, wp[i].error_estimate)
    for i in 1:3)
    @test all(startswith(plt[1][i + 3][:label], wp_names[i]) for i in 1:3)
    @test_throws BoundsError plt[1][7]
    @test_nowarn plot(wp, view = :dt_convergence, color = [:red :orange :green])
    @test_nowarn plot(
        wp, view = :dt_convergence, color = :lightblue, title = "Δt Convergence")
end

@testset "SDE WorkPrecisionSet" begin
    prob = remake(prob_sde_additivesystem, tspan = (0.0, 1.0))

    reltols = 1.0 ./ 10.0 .^ (1:5)
    abstols = reltols#[0.0 for i in eachindex(reltols)]
    setups = [Dict(:alg => SRIW1())
              Dict(:alg => EM(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1))
              Dict(:alg => RKMil(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => SRIW1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => SRA1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => SRA1())]
    names = ["SRIW1", "EM", "RKMil", "SRIW1 Fixed", "SRA1 Fixed", "SRA1"]
    test_dt = 0.1
    wp = WorkPrecisionSet(prob, abstols, reltols, setups; numruns = 10,
        names = names, maxiters = 1e7, error_estimate = :l2)

    plt = @test_nowarn plot(wp)
    for i in 1:length(names)
        @test plt[1][i][:x] ≈ getproperty(wp[i].errors, wp[i].error_estimate)
        @test plt[1][i][:label] == names[i]
    end
end

@testset "Analyticless SDE WorkPrecisionSet" begin
    prob = remake(prob_sde_additivesystem, tspan = (0.0, 1.0))

    reltols = 1.0 ./ 10.0 .^ (1:5)
    abstols = reltols#[0.0 for i in eachindex(reltols)]
    setups = [Dict(:alg => SRIW1())
              Dict(:alg => EM(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1))
              Dict(:alg => RKMil(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => SRIW1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => SRA1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => SRA1())]
    names = ["SRIW1", "EM", "RKMil", "SRIW1 Fixed", "SRA1 Fixed", "SRA1"]
    test_dt = 0.1
    wp = WorkPrecisionSet(prob, abstols, reltols, setups, test_dt;
        numruns = 5, names = names, error_estimate = :l2,
        appxsol_setup = Dict(:alg => RKMilGeneral(; ii_approx = IICommutative())), maxiters = 1e7)

    plt = @test_nowarn plot(wp)
    for i in 1:length(names)
        @test plt[1][i][:x] ≈ getproperty(wp[i].errors, wp[i].error_estimate)
        @test plt[1][i][:label] == names[i]
    end
end

@testset failfast=true "Complex SDE WorkPrecisionSet" begin
    # Linear SDE system
    f_lin = function (du, u, p, t)
        du = -0.5 .* u
    end

    g_lin = function (du, u, p, t)
        du = im .* u
    end

    lin_analytic = function (u₀, p, t, Wt)
        u₀ .* exp.(im .* Wt)
    end

    tspan = (0.0, 1.0)
    noise = StochasticDiffEq.RealWienerProcess!(0.0, [0.0], [0.0])
    prob = SDEProblem(SDEFunction(f_lin, g_lin; analytic = lin_analytic),
        ComplexF64[1.0], tspan, noise = noise)

    reltols = 1.0 ./ 10.0 .^ (1:5)
    abstols = reltols#[0.0 for i in eachindex(reltols)]
    setups = [Dict(:alg => EM(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1))
              Dict(:alg => RKMil(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => SRA1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => SRA1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => true)]
    names = ["EM", "RKMil", "SRA1 Fixed", "SRA1"]
    wp = WorkPrecisionSet(prob, abstols, reltols, setups; numruns = 10, numruns_error = 1,
        names = names, maxiters = 1e7, error_estimate = :l2)

    plt = @test_nowarn plot(wp)
    for i in 1:length(names)
        @test plt[1][i][:x] ≈ getproperty(wp[i].errors, wp[i].error_estimate)
        @test plt[1][i][:label] == names[i]
    end
end

@testset failfast=true "Non-diagonal SDE WorkPrecisionSet" begin
    # Linear SDE system
    f_lin = function (du, u, p, t)
        du = -0.5 .* u
    end

    g_lin = function (du, u, p, t)
        du[1, 1] = im
        du[2, 1] = im
        du[3, 1] = 0.1
        du[1, 2] = 0.1
        du[2, 2] = 0.1
        du[3, 2] = 0.2
    end

    tspan = (0.0, 1.0)
    noise_rate_prototype = zeros(ComplexF64, 3, 2)
    noise = StochasticDiffEq.RealWienerProcess!(0.0, [0.0, 0.0], [0.0, 0.0])
    prob = SDEProblem(SDEFunction(f_lin, g_lin),
        ComplexF64[1.0, 0.0, 0.0], tspan, noise = noise, noise_rate_prototype = noise_rate_prototype)

    reltols = 1.0 ./ 10.0 .^ (1:5)
    abstols = reltols#[0.0 for i in eachindex(reltols)]
    setups = [Dict(:alg => EM(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1))
              Dict(:alg => RKMilGeneral(; ii_approx = IICommutative()),
                  :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false);
              Dict(:alg => EulerHeun(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
                  :adaptive => false)
              Dict(:alg => LambaEulerHeun(),
                  :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1), :adaptive => true)]
    names = ["EM", "RKMilGeneral", "EulerHeun Fixed", "LambaEulerHeun"]
    test_dt = 0.1#(1.0 / 5.0)^6
    wp = WorkPrecisionSet(prob, abstols, reltols, setups, test_dt;
        numruns = 5, names = names, error_estimate = :l2,
        appxsol_setup = Dict(:alg => RKMilGeneral(; ii_approx = IICommutative())), maxiters = 1e7)

    plt = @test_nowarn plot(wp)
    for i in 1:length(names)
        @test plt[1][i][:x] ≈ getproperty(wp[i].errors, wp[i].error_estimate)
        @test plt[1][i][:label] == names[i]
    end
end