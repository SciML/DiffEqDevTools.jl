using DiffEqDevTools
using Base.Test

# write your own tests here
println("Benchmark Tests")
@time @test include("benchmark_tests.jl")
