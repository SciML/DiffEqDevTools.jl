"""
TestSolution

"""
type TestSolution <: DESolution
  u
  interp
  dense
end
(T::TestSolution)(t) = T.interp(t)
TestSolution(u) = TestSolution(u,nothing,false)
TestSolution(u,interp) = TestSolution(u,interp,true)
TestSolution(interp::DESolution) = TestSolution(nothing,interp,true)

"""
`appxtrue!(sol::AbstractODESolution,sol2::TestSolution)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.
"""
function appxtrue!(sol::AbstractODESolution,sol2::TestSolution)
  if sol2.u == nothing && sol2.dense
    sol2.u = sol2(sol.t[end])
  end
  errors = Dict(:final=>mean(abs.(sol[end]-sol2[end])))
  if sol2.dense
    timeseries_analytic = sol2(sol.t)
    errors = Dict(:final=>mean(abs.(sol[end]-sol2[end])),:l∞=>maximum(vecvecapply((x)->abs.(x),sol[:]-timeseries_analytic)),:l2=>sqrt(mean(vecvecapply((x)->float(x).^2,sol[:]-timeseries_analytic))))
    if !(typeof(sol) <: SDESolution) && sol.dense
      densetimes = collect(linspace(sol.t[1],sol.t[end],100))
      interp_u = sol(densetimes)
      interp_analytic = sol2(densetimes)
      interp_errors = Dict(:L∞=>maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic)),:L2=>sqrt(mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic))))
      errors = merge(errors,interp_errors)
    end
  end
  sol.u_analytic = sol2.u
  sol.errors = errors
  nothing
end

"""
`appxtrue!(sol::FEMSolution,sol2::FEMSolution)`

Adds the solution from `sol2` to the `FEMSolution` object `sol`.
Useful to add a quasi-true solution when none is known by
computing once at a very small time/space step and taking
that solution as the "true" solution
"""
function appxtrue!(sol::AbstractFEMSolution,sol2::AbstractFEMSolution)
  sol.u_analytic = sol2.u
  sol.errors = Dict(:l∞=>maximum(abs.(sol.u-sol.u_analytic)),:l2=>norm(sol.u-sol.u_analytic,2))
  nothing
end

"""
`appxtrue!(sol::AbstractODESolution,sol2::AbstractODESolution)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.
"""
function appxtrue!(sol::AbstractODESolution,sol2::AbstractODESolution)
  errors = Dict(:final=>mean(abs.(sol[end]-sol2[end])))
  if !(typeof(sol2) <: SDESolution) && sol2.dense
    timeseries_analytic = sol2(sol.t)
    errors = Dict(:final=>mean(abs.(sol[end]-sol2[end])),:l∞=>maximum(vecvecapply((x)->abs.(x),sol[:]-timeseries_analytic)),:l2=>sqrt(mean(vecvecapply((x)->float(x).^2,sol[:]-timeseries_analytic))))
    if !(typeof(sol) <: SDESolution) && sol.dense
      densetimes = collect(linspace(sol.t[1],sol.t[end],100))
      interp_u = sol(densetimes)
      interp_analytic = sol2(densetimes)
      interp_errors = Dict(:L∞=>maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic)),:L2=>sqrt(mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic))))
      errors = merge(errors,interp_errors)
    end
  end

  sol.u_analytic = sol2.u
  sol.errors = errors
  nothing
end
