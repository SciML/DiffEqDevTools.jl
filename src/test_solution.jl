"""
TestSolution

"""
mutable struct TestSolution{T, N, hasinterp, tType, uType, iType} <:
    AbstractTimeseriesSolution{T, N, uType}
    t::tType
    u::uType
    interp::iType
    dense::Bool
    retcode::ReturnCode.T
end
(T::TestSolution)(t) = T.interp(t)
function TestSolution(t, u)
    T = eltype(eltype(u))
    N = length((size(u[1])..., length(u)))
    return TestSolution{T, N, false, typeof(t), typeof(u), Nothing}(
        t, u, nothing, false,
        ReturnCode.Success
    )
end
function TestSolution(t, u, interp)
    T = eltype(eltype(u))
    N = length((size(u[1])..., length(u)))
    return TestSolution{T, N, true, typeof(t), typeof(u), typeof(interp)}(
        t, u, interp, true,
        ReturnCode.Success
    )
end
function TestSolution(interp::DESolution)
    return TestSolution{Nothing, 0, true, Nothing, Nothing, typeof(interp)}(
        nothing, nothing,
        interp, true,
        ReturnCode.Success
    )
end
function hasinterp(
        ::TestSolution{
            T, N, hi, tType, uType, iType,
        }
    ) where {
        T, N, hi, tType,
        uType, iType,
    }
    return hi
end
"""
    default_reduction(x)

Default reduction function for error calculation. Returns the input unchanged,
allowing the standard error metrics to compute errors over all components.
"""
default_reduction(x) = x

"""
`appxtrue(sol::AbstractODESolution,sol2::TestSolution; reduction=default_reduction)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.

## Arguments
- `sol`: The solution to compute errors for
- `sol2`: The reference (higher-order) solution
- `reduction`: A function to reduce/transform the solution difference before
  computing error metrics. Defaults to `default_reduction` (identity).
  For example, to compute error only on the third component of a system,
  use `reduction = x -> x[3]`. To compute error on a function of the solution,
  use `reduction = x -> some_function(x)`.
"""
function appxtrue(sol::AbstractODESolution, sol2::TestSolution; reduction = default_reduction)
    if sol2.u == nothing && hasinterp(sol2)
        _sol = TestSolution(sol.t, sol2(sol.t).u, sol2)
    else
        _sol = sol2
    end

    # Apply reduction to the final difference
    final_diff = reduction(sol.u[end] - _sol.u[end])
    errors = Dict(:final => recursive_mean(abs.(final_diff)))
    if _sol.dense
        timeseries_analytic = _sol(sol.t)
        # Apply reduction to each difference in the timeseries
        reduced_diff = vecvecapply(reduction, sol - timeseries_analytic)
        errors[:l∞] = maximum(vecvecapply((x) -> abs.(x), reduced_diff))
        errors[:l2] = sqrt(
            recursive_mean(
                vecvecapply(
                    (x) -> float(x) .^ 2,
                    reduced_diff
                )
            )
        )
        densetimes = collect(range(sol.t[1], stop = sol.t[end], length = 100))
        interp_u = sol(densetimes)
        interp_analytic = _sol(densetimes)
        # Apply reduction to dense interpolation differences
        reduced_interp_diff = vecvecapply(reduction, interp_u - interp_analytic)
        interp_errors = Dict(
            :L∞ => maximum(
                vecvecapply(
                    (x) -> abs.(x),
                    reduced_interp_diff
                )
            ),
            :L2 => sqrt(
                recursive_mean(
                    vecvecapply(
                        (x) -> float(x) .^ 2,
                        reduced_interp_diff
                    )
                )
            )
        )
        errors = merge(errors, interp_errors)
    else
        timeseries_analytic = sol2.u
        if sol.t == sol2.t
            # Apply reduction to each difference in the timeseries
            reduced_diff = vecvecapply(reduction, sol - timeseries_analytic)
            errors[:l∞] = maximum(vecvecapply((x) -> abs.(x), reduced_diff))
            errors[:l2] = sqrt(
                recursive_mean(
                    vecvecapply(
                        (x) -> float(x) .^ 2,
                        reduced_diff
                    )
                )
            )
        end
    end
    return DiffEqBase.build_solution(sol, timeseries_analytic, errors)
end

"""
`appxtrue(sol::AbstractODESolution,sol2::AbstractODESolution; reduction=default_reduction)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.

## Arguments
- `sol`: The solution to compute errors for
- `sol2`: The reference (higher-order) solution
- `timeseries_errors`: Whether to compute timeseries errors (default: sol2.dense)
- `dense_errors`: Whether to compute dense interpolation errors (default: sol2.dense)
- `reduction`: A function to reduce/transform the solution difference before
  computing error metrics. Defaults to `default_reduction` (identity).
  For example, to compute error only on the third component of a system,
  use `reduction = x -> x[3]`. To compute error on a function of the solution,
  use `reduction = x -> some_function(x)`.
"""
function appxtrue(
        sol::AbstractODESolution, sol2::AbstractODESolution;
        timeseries_errors::Bool = sol2.dense, dense_errors::Bool = sol2.dense,
        reduction = default_reduction
    )
    # Apply reduction to the final difference
    final_diff = reduction(sol.u[end] - sol2.u[end])
    errors = Dict(:final => recursive_mean(abs.(final_diff)))
    if sol2.dense
        timeseries_analytic = sol2(sol.t)
        # Apply reduction to each difference in the timeseries
        reduced_diff = vecvecapply(reduction, sol - timeseries_analytic)
        errors[:l∞] = maximum(vecvecapply((x) -> abs.(x), reduced_diff))
        errors[:l2] = sqrt(
            recursive_mean(
                vecvecapply(
                    (x) -> float(x) .^ 2,
                    reduced_diff
                )
            )
        )
        if dense_errors
            densetimes = collect(range(sol.t[1], stop = sol.t[end], length = 100))
            interp_u = sol(densetimes)
            interp_analytic = sol2(densetimes)
            # Apply reduction to dense interpolation differences
            reduced_interp_diff = vecvecapply(reduction, interp_u - interp_analytic)
            interp_errors = Dict(
                :L∞ => maximum(
                    vecvecapply(
                        (x) -> abs.(x),
                        reduced_interp_diff
                    )
                ),
                :L2 => sqrt(
                    recursive_mean(
                        vecvecapply(
                            (x) -> float(x) .^ 2,
                            reduced_interp_diff
                        )
                    )
                )
            )
            errors = merge(errors, interp_errors)
        end
    else
        timeseries_analytic = sol2.u
        if timeseries_errors && sol.t == sol2.t
            # Apply reduction to each difference in the timeseries
            reduced_diff = vecvecapply(reduction, sol - timeseries_analytic)
            errors[:l∞] = maximum(vecvecapply((x) -> abs.(x), reduced_diff))
            errors[:l2] = sqrt(
                recursive_mean(
                    vecvecapply(
                        (x) -> float(x) .^ 2,
                        reduced_diff
                    )
                )
            )
        end
    end
    return DiffEqBase.build_solution(sol, timeseries_analytic, errors)
end

function appxtrue(sim::EnsembleSolution, appx_setup; kwargs...)
    _new_sols = Vector{DESolution}(length(sim.u))
    for i in eachindex(sim)
        prob = sim[i].prob
        prob2 = SDEProblem(
            prob.f, prob.g, prob.u0, prob.tspan,
            noise = NoiseWrapper(sim.u[i].W)
        )
        true_sol = solve(prob2, appx_setup[:alg]; appx_setup...)
        _new_sols[i] = appxtrue(sim.u[i], true_sol)
    end
    new_sols = convert(Vector{typeof(_new_sols[1])}, _new_sols)
    DiffEqBase.calculate_ensemble_errors(
        new_sols; converged = sim.converged,
        elapsedTime = sim.elapsedTime, kwargs...
    )
end
