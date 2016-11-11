using DiffEqDevTools
using Base.Test

# write your own tests here
println("Benchmark Tests")
@time @test include("benchmark_tests.jl")
println("ODE AppxTrue Tests")
@time @test include("ode_appxtrue_tests.jl")
println("ODE Tableau Convergence Tests")
!is_windows() && @time @test include("ode_tableau_convergence_tests.jl") ## Windows 32-bit fails on Butcher62 convergence test
