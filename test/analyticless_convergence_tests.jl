using OrdinaryDiffEq, StochasticDelayDiffEq, ParameterizedFunctions, Test, Random
using ParameterizedFunctions.ModelingToolkit # macro hygiene
f = @ode_def LotkaVolterra begin
  dx = 1.5x - x*y
  dy = -3y + x*y
end

prob = ODEProblem(f,big.([1.0;1.0]),(big(0.0),big(10.0)))

using DiffEqDevTools

dts = big(1/2).^(8:-1:4)
test_setup = Dict(:alg=>Vern9(),:reltol=>1e-25,:abstol=>1e-25)
sim1 = analyticless_test_convergence(dts,prob,Tsit5(),test_setup)
sim2 = analyticless_test_convergence(dts,prob,Vern9(),test_setup)

@test abs(sim1.𝒪est[:final]-5) < 0.2
@test abs(sim2.𝒪est[:final]-9) < 0.2

function f2(du,u,p,t)
  du .= 1.01u
end
function g2(du,u,p,t)
  du .= 1.01u
end
prob = SDEProblem(f2,g2,[1.0;1.0],(0.0,1.0))

using StochasticDiffEq

dts = (1/2).^(7:-1:3)
test_dt = 1/2^8
Random.seed!(100)
sim1 = analyticless_test_convergence(dts,prob,SRIW1(),test_dt,trajectories=100)
@test abs(sim1.𝒪est[:final]-1.5) < 0.4
@show sim1.𝒪est[:final]

dts = (1/2).^(7:-1:4)
test_dt = 1/2^8
sim2 = analyticless_test_convergence(dts,prob,SRIW1(),test_dt,trajectories=100, use_noise_grid=false)
@test abs(sim2.𝒪est[:final]-1.5) < 0.3
@show sim2.𝒪est[:final]

### SDDE

function hayes_modelf(du,u,h,p,t)
    τ,a,b,c,α,β,γ = p
    du .= a.*u .+ b .* h(p,t-τ) .+ c
end
function hayes_modelg(du,u,h,p,t)
    τ,a,b,c,α,β,γ = p
    du .= α.*u .+ β.*h(p,t-τ) .+ γ
end
h(p,t) = (ones(1) .+ t);
tspan = (0.,10.)

pmul = [1.0,-4.,-2.,10.,-1.3,-1.2, 1.1]
padd = [1.0,-4.,-2.,10.,-0.0,-0.0, 0.1]

prob = SDDEProblem(hayes_modelf, hayes_modelg, [1.], h, tspan, pmul; constant_lags = (pmul[1],));
dts = (1/2).^(7:-1:3)
test_dt = 1/2^8
sim2 = analyticless_test_convergence(dts,prob,RKMil(),test_dt,trajectories=100, use_noise_grid=false)
@test abs(sim2.𝒪est[:final]-1.0) < 0.3
