using OrdinaryDiffEq, ParameterizedFunctions, Base.Test

f = @ode_def_nohes LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=1.5 b=1 c=3 d=1

prob = ODEProblem(f,big.([1.0;1.0]),(big(0.0),big(10.0)))

using DiffEqDevTools

dts = big(1./2.).^(8:-1:4)
test_setup = Dict(:alg=>Vern9(),:reltol=>1e-25,:abstol=>1e-25)
sim1 = analyticless_test_convergence(dts,prob,Tsit5(),test_setup)
sim2 = analyticless_test_convergence(dts,prob,Vern9(),test_setup)

@test sim1.𝒪est[:final]-5 < 0.2
@test sim2.𝒪est[:final]-9 < 0.2

function f2(t,u,du)
  du .= 1.01u
end
function g2(t,u,du)
  du .= 1.01u
end
prob = SDEProblem(f2,g2,[1.0;1.0],(0.0,10.0))

using StochasticDiffEq

dts = 1./2.^(6:-1:3)
test_dt = 1/2^8
sim1 = analyticless_test_convergence(dts,prob,SRIW1(),test_dt,numMonte=200)
@test sim1.𝒪est[:final]-1.5 < 0.3
