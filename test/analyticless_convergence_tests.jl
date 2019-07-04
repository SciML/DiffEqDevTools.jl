using OrdinaryDiffEq, ParameterizedFunctions, Test, Random

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

@test sim1.𝒪est[:final]-5 < 0.2
@test sim2.𝒪est[:final]-9 < 0.2

function f2(du,u,p,t)
  du .= 1.01u
end
function g2(du,u,p,t)
  du .= 1.01u
end
prob = SDEProblem(f2,g2,[1.0;1.0],(0.0,10.0))

using StochasticDiffEq

dts = (1/2).^(6:-1:3)
test_dt = 1/2^8
Random.seed!(100)
sim1 = analyticless_test_convergence(dts,prob,SRIW1(),test_dt,trajectories=400)
@test sim1.𝒪est[:final]-1.5 < 0.3
