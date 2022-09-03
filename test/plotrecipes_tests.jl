using Test
using OrdinaryDiffEq, DiffEqDevTools, Plots

using Random
Random.seed!(123)

gr()

# 2D Linear ODE
function f2d(du, u, p, t)
    @inbounds for i in eachindex(u)
        du[i] = cos(t) * u[i]
    end
end

function f2d_analytic(u₀, p, t)
    u₀ * exp(sin(t))
end

tspan = (0.0,10.0)
prob = ODEProblem(ODEFunction(f2d,analytic=f2d_analytic),rand(10,10),tspan)

abstols = 1.0 ./ 10.0 .^ (3:8)
reltols = 1.0 ./ 10.0 .^ (0:5)

setups = [
    Dict(:alg=>DP5())
    Dict(:alg=>Tsit5())
]
names = ["DP5", "Tsit5"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups,names=names,save_everystep=false,numruns=100)

plt = @test_nowarn plot(wp)
@test plt[1][1][:x] ≈ wp[1].errors
@test plt[1][2][:x] ≈ wp[2].errors
@test plt[1][1][:label] == names[1]
@test plt[1][2][:label] == names[2]
@test_throws ArgumentError plot(wp, view=:dt_convergence)

dts = 1.0./2.0.^((1:length(reltols)) .+ 1)
setups = [
    Dict(:alg=>Euler(),:dts=>dts)
    Dict(:alg=>Heun(),:dts=>dts)
    Dict(:alg=>DP5(),:dts=>dts,:adaptive=>false)
    Dict(:alg=>Tsit5(),:dts=>dts,:adaptive=>false)
    Dict(:alg=>Tsit5())
]
names = ["Euler", "Heun", "DP5 fixed step", "Tsit5 fixed step", "Tsit5 adaptive"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups,names=names,save_everystep=false,numruns=100)

plt = @test_nowarn plot(wp)
@test all(plt[1][i][:x] ≈ wp[i].errors for i in 1:5)
@test all(plt[1][i][:label] == names[i] for i in 1:5)

plt = @test_nowarn plot(wp, view=:dt_convergence, legend=:bottomright)
@test all(plt[1][i][:x] == plt[1][i+4][:x] == dts == wp.setups[i][:dts] for i in 1:4)
@test all(plt[1][i+4][:y] ≈ wp[i].errors for i in 1:4)
@test all(startswith(plt[1][i+4][:label], names[i]) for i in 1:4)
@test_throws BoundsError plt[1][9]
