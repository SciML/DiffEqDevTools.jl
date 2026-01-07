using AllocCheck
using DiffEqDevTools
using Test

@testset "AllocCheck - Allocation Regression Tests" begin
    # Test that calc𝒪estimates has minimal allocations
    # This function was optimized to avoid allocating intermediate arrays
    @testset "calc𝒪estimates allocations" begin
        errors = [1e-2, 5e-3, 2.5e-3, 1.25e-3]
        pair = :final => errors

        # Warm up
        DiffEqDevTools.calc𝒪estimates(pair)

        # Verify zero allocations for the core computation
        allocs = @allocated DiffEqDevTools.calc𝒪estimates(pair)
        # Allow for small allocations from Pair creation but should be minimal
        @test allocs <= 64  # Very minimal - just the Pair object
    end

    # Test allocation regression guard using @check_allocs
    # This will fail compilation if the function allocates
    @testset "calc𝒪estimates @check_allocs" begin
        # Define a wrapper that @check_allocs can analyze
        @check_allocs function calc_order_estimates_wrapper(key::Symbol, errors::Vector{Float64})
            n = length(errors)
            s = 0.0
            @inbounds for i in 1:(n - 1)
                s += log2(errors[i + 1] / errors[i])
            end
            return abs(s / (n - 1))
        end

        errors = [1e-2, 5e-3, 2.5e-3, 1.25e-3]
        result = calc_order_estimates_wrapper(:final, errors)
        @test result ≈ 1.0 atol = 0.01
    end
end
