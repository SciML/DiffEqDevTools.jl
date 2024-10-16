using Random
using StochasticDiffEq, DiffEqDevTools, Test
using SDEProblemLibrary: prob_sde_additivesystem

prob = prob_sde_additivesystem
prob = SDEProblem(prob.f, prob.g, prob.u0, (0.0, 0.1), prob.p)

reltols = 1.0 ./ 10.0 .^ (1:4)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [Dict(:alg => SRIW1())
          Dict(:alg => EM(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
              :adaptive => false)
          Dict(:alg => RKMil(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
              :adaptive => false)
          Dict(:alg => SRIW1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
              :adaptive => false)
          Dict(:alg => SRA1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
              :adaptive => false)
          Dict(:alg => SRA1())]
_names = ["SRIW1", "EM", "RKMil", "SRIW1 Fixed", "SRA1 Fixed", "SRA1"]
test_dt = 0.1
# wp = WorkPrecisionSet(prob, abstols, reltols, setups, test_dt;
#     numruns = 5, names = _names, error_estimate = :l2)

# se = get_sample_errors(prob, setups[1], numruns = 100, solution_runs = 100)
# se = get_sample_errors(prob, setups[1], numruns = [5, 10, 25, 50, 100, 1000],
#     solution_runs = 100)

println("Now weak error without analytical solution")

prob2 = SDEProblem((du, u, p, t) -> prob.f(du, u, p, t), prob.g, prob.u0, (0.0, 0.1),
    prob.p)
test_dt = 1 / 10^4
appxsol_setup = Dict(:alg => SRIW1(), :abstol => 1e-4, :reltol => 1e-4)

wp = WorkPrecisionSet(prob2, abstols, reltols, setups, test_dt;
    appxsol_setup = appxsol_setup,
    numruns = 5, names = _names, error_estimate = :l2)

using Test
using OrdinaryDiffEq, StochasticDiffEq, DiffEqDevTools, Plots
import SDEProblemLibrary: prob_sde_additivesystem

using Random
Random.seed!(123)

gr()

@testset "Analyticless SDE WorkPrecisionSet" begin
    prob0 = prob_sde_additivesystem
    prob = SDEProblem((du, u, p, t) -> prob0.f(du, u, p, t), prob0.g, prob0.u0, (0.0, 0.1),
        prob0.p)

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
        appxsol_setup = Dict(:alg => RKMilCommute(), :abstol => 1e-4, :reltol => 1e-4))

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