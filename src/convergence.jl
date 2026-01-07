"""
    ConvergenceSimulation{SolType}

Stores the results of a convergence test across multiple timesteps or parameter values.

# Fields
- `solutions`: Array of solutions at different parameter values
- `errors`: Dictionary of error estimates for different error types
- `N`: Number of simulations performed
- `auxdata`: Auxiliary data (e.g., timesteps used)
- `𝒪est`: Dictionary of estimated convergence orders for each error type
- `convergence_axis`: The parameter being varied (e.g., timesteps)

# Indexing
Supports array-like indexing to access individual solutions:
- `sim[i]` returns the i-th solution
- `sim[i, j, ...]` passes additional indices to the i-th solution
"""
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
    #𝒪est = Dict(map(calc𝒪estimates,errors))
    𝒪esttmp = Dict() #Makes Dict of Any to be more compatible
    for (k, v) in 𝒪est
        if length(v) == 1
            push!(𝒪esttmp, Pair(k, v[1]))
        else
            push!(𝒪esttmp, Pair(k, v))
        end
    end
    𝒪est = 𝒪esttmp
    return (ConvergenceSimulation(solutions, errors, N, auxdata, 𝒪est, convergence_axis))
end

"""
    test_convergence(dts, prob::Union{AbstractRODEProblem, AbstractSDEProblem, AbstractEnsembleProblem},
                     alg, ensemblealg=EnsembleThreads(); trajectories, kwargs...)

Test the convergence rate of a stochastic differential equation solver by running
Monte Carlo simulations at different timesteps.

# Arguments
- `dts`: Array of timesteps to test
- `prob`: The SDE/RODE problem or EnsembleProblem to solve
- `alg`: The algorithm to test
- `ensemblealg`: Parallelization strategy (default: `EnsembleThreads()`)

# Keyword Arguments
- `trajectories`: Number of Monte Carlo trajectories to run (required)
- `save_start=true`: Save the initial condition
- `save_everystep=true`: Save at every timestep
- `timeseries_errors=save_everystep`: Calculate timeseries error metrics
- `adaptive=false`: Use adaptive timestepping
- `weak_timeseries_errors=false`: Calculate weak error metrics over time
- `weak_dense_errors=false`: Calculate weak error metrics on dense output
- `expected_value=nothing`: Expected value for weak error calculation (if known analytically)

# Returns
A `ConvergenceSimulation` object containing:
- Solutions at each timestep
- Error estimates (`:l2`, `:l∞`, `:L2`, `:L∞`, `:weak_final`, etc.)
- Estimated convergence orders in the `𝒪est` field

# Example
```julia
using StochasticDiffEq, DiffEqDevTools

f(du, u, p, t) = (du .= 1.01u)
g(du, u, p, t) = (du .= 0.87u)
prob = SDEProblem(f, g, [1.0], (0.0, 1.0))

dts = (1 / 2) .^ (5:-1:3)
sim = test_convergence(dts, prob, SRIW1(), trajectories = 1000)
println("Estimated order: ", sim.𝒪est[:final])
```
"""
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
    analyticless_test_convergence(dts, prob::Union{AbstractRODEProblem, AbstractSDEProblem, AbstractSDDEProblem},
                                   alg, test_dt; trajectories=100, kwargs...)

Test convergence of a stochastic solver without an analytical solution by using
a high-resolution reference solution computed at `test_dt`.

# Arguments
- `dts`: Array of timesteps to test
- `prob`: The SDE/RODE/SDDE problem to solve
- `alg`: The algorithm to test
- `test_dt`: Small timestep for computing reference solutions

# Keyword Arguments
- `trajectories=100`: Number of Monte Carlo trajectories
- `save_everystep=true`: Save at every timestep
- `timeseries_errors=save_everystep`: Calculate timeseries errors
- `adaptive=false`: Use adaptive timestepping
- `weak_timeseries_errors=false`: Calculate weak timeseries errors
- `weak_dense_errors=false`: Calculate weak dense errors
- `use_noise_grid=true`: Use noise grid for reproducibility
- `verbose=true`: Print progress information

# Returns
A `ConvergenceSimulation` object with error estimates and convergence orders.

# Example
```julia
using StochasticDiffEq, DiffEqDevTools

f(du, u, p, t) = (du .= 1.01u)
g(du, u, p, t) = (du .= 0.87u)
prob = SDEProblem(f, g, [1.0], (0.0, 1.0))

dts = (1 / 2) .^ (7:-1:4)
test_dt = 1 / 2^8  # Fine timestep for reference
sim = analyticless_test_convergence(dts, prob, SRIW1(), test_dt, trajectories = 100)
```

# Notes
This function generates a reference solution at `test_dt` for each trajectory,
then compares solutions at coarser timesteps to estimate convergence rates.
"""
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
    for j in 1:trajectories
        if verbose
            @info "Monte Carlo iteration: $j/$trajectories"
        end
        t = prob.tspan[1]:test_dt:prob.tspan[2]
        if use_noise_grid
            T = eltype(prob.u0)
            if prob.noise_rate_prototype === nothing
                brownian_values = cumsum(
                    [
                        [zeros(T, size(prob.u0))];
                        [
                            sqrt(test_dt) * randn(T, size(prob.u0))
                                for i in 1:(length(t) - 1)
                        ]
                    ]
                )
                brownian_values2 = cumsum(
                    [
                        [zeros(T, size(prob.u0))];
                        [
                            sqrt(test_dt) * randn(T, size(prob.u0))
                                for i in 1:(length(t) - 1)
                        ]
                    ]
                )
            else
                brownian_values = cumsum(
                    [
                        [zeros(T, size(prob.noise_rate_prototype, 2))];
                        [
                            sqrt(test_dt) *
                                randn(T, size(prob.noise_rate_prototype, 2))
                                for i in 1:(length(t) - 1)
                        ]
                    ]
                )
                brownian_values2 = cumsum(
                    [
                        [zeros(T, size(prob.noise_rate_prototype, 2))];
                        [
                            sqrt(test_dt) *
                                randn(T, size(prob.noise_rate_prototype, 2))
                                for i in 1:(length(t) - 1)
                        ]
                    ]
                )
            end
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

"""
    test_convergence(dts, prob::Union{AbstractODEProblem, AbstractDAEProblem}, alg; kwargs...)

Test the convergence rate of an ODE/DAE solver at different timesteps. Requires
the problem to have an analytical solution defined.

# Arguments
- `dts`: Array of timesteps to test
- `prob`: The ODE or DAE problem with analytical solution
- `alg`: The algorithm to test

# Keyword Arguments
- `save_everystep=true`: Save the solution at every timestep
- `adaptive=false`: Use adaptive timestepping (typically false for convergence tests)

# Returns
A `ConvergenceSimulation` object containing solutions and convergence estimates.

# Example
```julia
using OrdinaryDiffEq, DiffEqDevTools

# Problem with analytical solution
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)

dts = 1 ./ 2 .^ (6:10)
sim = test_convergence(dts, prob, Tsit5())
println("Estimated order: ", sim.𝒪est[:final])  # Should be ≈5
```

# Notes
If the problem doesn't have an analytical solution, use `analyticless_test_convergence` instead.
"""
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

"""
    analyticless_test_convergence(dts, prob::AbstractODEProblem, alg, appxsol_setup; kwargs...)

Test convergence of an ODE solver without an analytical solution by comparing
against a high-accuracy reference solution.

# Arguments
- `dts`: Array of timesteps to test
- `prob`: The ODE problem (without analytical solution)
- `alg`: The algorithm to test
- `appxsol_setup`: Dictionary specifying the reference solution algorithm and tolerances,
  e.g., `Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)`

# Keyword Arguments
- `save_everystep=true`: Save at every timestep
- `adaptive=false`: Use adaptive timestepping (typically false for convergence tests)

# Returns
A `ConvergenceSimulation` object with error estimates and convergence orders.

# Example
```julia
using OrdinaryDiffEq, DiffEqDevTools

function lotka_volterra(du, u, p, t)
    du[1] = 1.5 * u[1] - u[1] * u[2]
    du[2] = -3 * u[2] + u[1] * u[2]
end
prob = ODEProblem(lotka_volterra, [1.0, 1.0], (0.0, 10.0))

dts = 1 ./ 2 .^ (6:9)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, Tsit5(), test_setup)
println("Estimated order: ", sim.𝒪est[:final])
```
"""
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
    error = error.second

    if ndims(error) > 1
        error = mean(error, 1)
    end
    S = Vector{eltype(error)}(undef, length(error) - 1)
    for i in 1:(length(error) - 1)
        S[i] = log2(error[i + 1] / error[i])
    end
    return (Pair(key, abs.(mean(S, dims = 1))))
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
