mutable struct ConvergenceSimulation{SolType}
    solutions::Array{SolType}
    errors::Any
    N::Any
    auxdata::Any
    𝒪est::Any
    convergence_axis::Any
end

function ConvergenceSimulation(
        solutions, convergence_axis;
        auxdata = nothing, additional_errors = nothing,
        expected_value = nothing
    )
    N = size(solutions, 1)
    uEltype = eltype(solutions[1].u[1])
    errors = Dict() #Should add type information
    if expected_value == nothing
        if isnothing(solutions[1].errors) || isempty(solutions[1].errors)
            error("Errors dictionary is empty. No analytical solution set. If you used `test_convergence` you may consider `analyticless_test_convergence` instead.")
        end
        for k in keys(solutions[1].errors)
            errors[k] = [mean(sol.errors[k]) for sol in solutions]
        end
    end
    if additional_errors != nothing
        for k in keys(additional_errors)
            errors[k] = additional_errors[k]
        end
    end

    𝒪est = Dict((calc𝒪estimates(p) for p in pairs(errors)))
    return (ConvergenceSimulation(solutions, errors, N, auxdata, 𝒪est, convergence_axis))
end

function test_convergence(
        dts::AbstractArray,
        prob::Union{
            AbstractRODEProblem, AbstractSDEProblem,
            AbstractEnsembleProblem,
        },
        alg, ensemblealg = EnsembleThreads();
        trajectories, save_start = true, save_everystep = true,
        timeseries_errors = save_everystep, adaptive = false,
        weak_timeseries_errors = false, weak_dense_errors = false,
        expected_value = nothing, kwargs...
    )
    N = length(dts)

    if prob isa AbstractEnsembleProblem
        ensemble_prob = prob
    else
        ensemble_prob = EnsembleProblem(prob)
    end

    _solutions = Array{Any}(undef, length(dts))
    for i in 1:length(dts)
        sol = solve(
            ensemble_prob, alg, ensemblealg; dt = dts[i], adaptive = adaptive,
            save_start = save_start, save_everystep = save_everystep,
            timeseries_errors = timeseries_errors,
            weak_timeseries_errors = weak_timeseries_errors,
            weak_dense_errors = weak_dense_errors, trajectories = Int(trajectories),
            kwargs...
        )
        @info "dt: $(dts[i]) ($i/$N)"
        _solutions[i] = sol
    end

    auxdata = Dict("dts" => dts)

    if expected_value == nothing
        solutions = [
            DiffEqBase.calculate_ensemble_errors(
                    sim;
                    weak_timeseries_errors = weak_timeseries_errors,
                    weak_dense_errors = weak_dense_errors
                )
                for sim in _solutions
        ]
        # Now Calculate Weak Errors
        additional_errors = Dict()
        for k in keys(solutions[1].weak_errors)
            additional_errors[k] = [sol.weak_errors[k] for sol in solutions]
        end

    else
        additional_errors = Dict()
        if length(expected_value) == 1 || expected_value isa Number
            additional_errors[:weak_final] = []
        else
            additional_errors[:weak_l2] = []
        end
        for sol in _solutions
            if length(expected_value) == 1 || expected_value isa Number
                weak_final = LinearAlgebra.norm(Statistics.mean(sol.u .- expected_value))
                push!(additional_errors[:weak_final], weak_final)
            else
                weak_l2 = LinearAlgebra.norm(Statistics.mean(sol .- expected_value))
                push!(additional_errors[:weak_l2], weak_l2)
            end
        end
        solutions = _solutions
    end

    return ConvergenceSimulation(
        solutions, dts, auxdata = auxdata,
        additional_errors = additional_errors,
        expected_value = expected_value
    )
end

"""
    _generate_brownian_values!(brownian_values, test_dt, T, dims)

Pre-allocate and fill Brownian motion values in-place.
"""
function _generate_brownian_values!(brownian_values, test_dt, T, dims)
    sqrt_dt = sqrt(test_dt)
    len = length(brownian_values)
    @inbounds brownian_values[1] = zeros(T, dims)
    @inbounds for i in 2:len
        brownian_values[i] = brownian_values[i - 1] + sqrt_dt * randn(T, dims)
    end
    return brownian_values
end

function analyticless_test_convergence(
        dts::AbstractArray,
        prob::Union{
            AbstractRODEProblem, AbstractSDEProblem,
            AbstractSDDEProblem,
        },
        alg, test_dt; trajectories = 100,
        save_everystep = true,
        timeseries_errors = save_everystep, adaptive = false,
        weak_timeseries_errors = false,
        weak_dense_errors = false, use_noise_grid = true,
        verbose = true,
        kwargs...
    )
    _solutions = []
    tmp_solutions = Array{Any}(undef, trajectories, length(dts))

    # Pre-compute time grid (constant across all trajectories)
    t = prob.tspan[1]:test_dt:prob.tspan[2]
    len_t = length(t)
    T = eltype(prob.u0)

    # Determine Brownian motion dimensions
    if prob.noise_rate_prototype === nothing
        noise_dims = size(prob.u0)
    else
        noise_dims = (size(prob.noise_rate_prototype, 2),)
    end

    for j in 1:trajectories
        if verbose
            @info "Monte Carlo iteration: $j/$trajectories"
        end
        if use_noise_grid
            # Pre-allocate Brownian motion arrays
            brownian_values = Vector{Array{T, length(noise_dims)}}(undef, len_t)
            brownian_values2 = Vector{Array{T, length(noise_dims)}}(undef, len_t)

            # Generate Brownian motions in-place
            _generate_brownian_values!(brownian_values, test_dt, T, noise_dims)
            _generate_brownian_values!(brownian_values2, test_dt, T, noise_dims)

            np = NoiseGrid(t, brownian_values, brownian_values2)

            if prob isa AbstractSDDEProblem
                _prob = SDDEProblem(
                    prob.f, prob.g, prob.u0, prob.h, prob.tspan, prob.p,
                    noise = np,
                    noise_rate_prototype = prob.noise_rate_prototype,
                    constant_lags = prob.constant_lags,
                    dependent_lags = prob.dependent_lags,
                    neutral = prob.neutral,
                    order_discontinuity_t0 = prob.order_discontinuity_t0,
                    prob.kwargs...
                )
            else
                _prob = SDEProblem(
                    prob.f, prob.g, prob.u0, prob.tspan, prob.p,
                    noise = np,
                    noise_rate_prototype = prob.noise_rate_prototype
                )
            end

            true_sol = solve(_prob, alg; adaptive = adaptive, dt = test_dt)

            for i in 1:length(dts)
                sol = solve(_prob, alg; dt = dts[i], adaptive = adaptive)
                err_sol = appxtrue(sol, true_sol)
                tmp_solutions[j, i] = err_sol
            end
        else
            # using NoiseWrapper doesn't lead to constant true_sol
            true_sol = solve(
                prob, alg; adaptive = adaptive, dt = test_dt,
                save_noise = true
            )
            _sol = deepcopy(true_sol)
            W1 = NoiseWrapper(_sol.W)

            _prob = remake(
                prob, u0 = prob.u0, p = prob.p, tspan = prob.tspan, noise = W1,
                noise_rate_prototype = prob.noise_rate_prototype
            )

            for i in 1:length(dts)
                W1 = NoiseWrapper(_sol.W)
                _prob = remake(
                    prob, u0 = prob.u0, p = prob.p, tspan = prob.tspan,
                    noise = W1, noise_rate_prototype = prob.noise_rate_prototype
                )
                sol = solve(_prob, alg; dt = dts[i], adaptive = adaptive)
                err_sol = appxtrue(sol, true_sol)
                tmp_solutions[j, i] = err_sol
            end
        end
    end
    _solutions = [EnsembleSolution(tmp_solutions[:, i], 0.0, true) for i in 1:length(dts)]
    solutions = [
        DiffEqBase.calculate_ensemble_errors(
                sim;
                weak_timeseries_errors = weak_timeseries_errors,
                weak_dense_errors = weak_dense_errors
            )
            for sim in _solutions
    ]
    auxdata = Dict("dts" => dts)
    # Now Calculate Weak Errors
    additional_errors = Dict()
    for k in keys(solutions[1].weak_errors)
        additional_errors[k] = [sol.weak_errors[k] for sol in solutions]
    end
    return ConvergenceSimulation(
        solutions, dts, auxdata = auxdata,
        additional_errors = additional_errors
    )
end

function test_convergence(
        dts::AbstractArray,
        prob::Union{AbstractODEProblem, AbstractDAEProblem}, alg;
        save_everystep = true, adaptive = false, kwargs...
    )
    N = length(dts)
    solutions = [
        solve(
                prob, alg; dt = dts[i], save_everystep = save_everystep,
                adaptive = adaptive, kwargs...
            ) for i in 1:N
    ]
    auxdata = Dict(:dts => dts)
    return ConvergenceSimulation(solutions, dts, auxdata = auxdata)
end

function analyticless_test_convergence(
        dts::AbstractArray, prob::AbstractODEProblem,
        alg, appxsol_setup;
        save_everystep = true, adaptive = false, kwargs...
    )
    true_sol = solve(prob, appxsol_setup[:alg]; appxsol_setup...)
    N = length(dts)
    _solutions = [
        solve(
                prob, alg; dt = dts[i], save_everystep = save_everystep,
                adaptive = adaptive, kwargs...
            ) for i in 1:N
    ]
    solutions = [appxtrue(sol, true_sol) for sol in _solutions]
    auxdata = Dict(:dts => dts)
    return ConvergenceSimulation(solutions, dts, auxdata = auxdata)
end

function test_convergence(probs, convergence_axis, alg; kwargs...)
    return ConvergenceSimulation([solve(prob, alg; kwargs...) for prob in probs], convergence_axis)
end

function test_convergence(c::ConvergenceSetup, alg::DEAlgorithm; kwargs...)
    return test_convergence(c.probs, c.convergence_axis, alg; kwargs...)
end

function calc𝒪estimates(error::Pair)
    key = error.first
    error_vals = error.second

    if ndims(error_vals) > 1
        error_vals = mean(error_vals, dims = 1)
    end
    n = length(error_vals)
    # Compute mean of log2 ratios directly without allocating intermediate vector
    # Use Float64 accumulator to handle Rational inputs
    s = 0.0
    @inbounds for i in 1:(n - 1)
        s += log2(error_vals[i + 1] / error_vals[i])
    end
    return Pair(key, abs(s / (n - 1)))
end

"""
length(simres::ConvergenceSimulation)

Returns the number of simulations in the Convergence Simulation
"""
Base.length(sim::ConvergenceSimulation) = sim.N
Base.getindex(sim::ConvergenceSimulation, i::Int) = sim.solutions[i]
Base.getindex(sim::ConvergenceSimulation, i::Int, I::Int...) = sim.solutions[i][I]
Base.lastindex(sim::ConvergenceSimulation) = lastindex(sim.solutions)
Base.firstindex(sim::ConvergenceSimulation) = firstindex(sim.solutions)
Base.length(sim::ConvergenceSetup) = sim.probs
