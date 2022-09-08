using Test
using OrdinaryDiffEq, DiffEqDevTools, Plots

using Random
Random.seed!(123)

gr()

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
prob = ODEProblem(ODEFunction(f_rode_lin, analytic = f_rode_lin_analytic), rand(10, 10),
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
@test plt[1][1][:x] ≈ wp[1].errors
@test plt[1][2][:x] ≈ wp[2].errors
@test plt[1][1][:label] == wp_names[1]
@test plt[1][2][:label] == wp_names[2]
@test_nowarn plot(wp, color = [1 2])
@test_nowarn plot(wp, color = :blue)
@test_throws ArgumentError plot(wp, view = :dt_convergence)

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
@test all(plt[1][i][:x] ≈ wp[i].errors for i in 1:4)
@test all(plt[1][i][:label] == wp_names[i] for i in 1:4)

plt = @test_nowarn plot(wp, view = :dt_convergence, legend = :bottomright)
@test all(plt[1][i][:x] == plt[1][i + 3][:x] == dts == wp.setups[i][:dts] for i in 1:3)
@test all(plt[1][i + 3][:y] ≈ wp[i].errors for i in 1:3)
@test all(startswith(plt[1][i + 3][:label], wp_names[i]) for i in 1:3)
@test_throws BoundsError plt[1][7]
@test_nowarn plot(wp, view = :dt_convergence, color = [:red :orange :green])
@test_nowarn plot(wp, view = :dt_convergence, color = :lightblue, title = "Δt Convergence")
