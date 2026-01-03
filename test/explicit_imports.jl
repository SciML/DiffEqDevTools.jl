using ExplicitImports
using DiffEqDevTools
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(DiffEqDevTools) === nothing
    @test check_no_stale_explicit_imports(DiffEqDevTools) === nothing
end
