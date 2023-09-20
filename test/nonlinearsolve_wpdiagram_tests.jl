# Fetch pakages.
using NonlinearSolve, DiffEqDevTools, Plots

# Prepares NonlinearProblem.
let
    f(u, p) = u .* u .- p
    u0 = [1.0, 1.0]
    p = 2.0
    static_prob = NonlinearProblem(f, u0, p)
    real_sol = solve(static_prob, NewtonRaphson(), reltol = 1e-15, abstol = 1e-15)

    # Sets WP input.
    abstols = 1.0 ./ 10.0 .^ (8:12)
    reltols = 1.0 ./ 10.0 .^ (8:12)
    setups = [Dict(:alg => NewtonRaphson())
        Dict(:alg => TrustRegion())]
    solnames = ["NewtonRaphson"; "TrustRegion"]

    # Makes WP-diagram
    wp = WorkPrecisionSet(static_prob,
        abstols,
        reltols,
        setups;
        names = solnames,
        numruns = 100,
        appxsol = real_sol,
        error_estimate = :l2)

    # Checks that all errors are small (they definitely should be).
    all(vcat(getfield.(wp.wps, :errors)...) .< 10e-9)
end
