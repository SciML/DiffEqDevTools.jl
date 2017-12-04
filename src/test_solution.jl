"""
TestSolution

"""
type TestSolution{T,N,hasinterp,tType,uType,iType} <: AbstractTimeseriesSolution{T,N}
  t::tType
  u::uType
  interp::iType
  dense::Bool
  retcode::Symbol
end
(T::TestSolution)(t) = T.interp(t)
function TestSolution(t,u)
  T = eltype(eltype(u))
  N = length((size(u[1])..., length(u)))
  TestSolution{T,N,false,typeof(t),typeof(u),Void}(t,u,nothing,false,:Success)
end
function TestSolution(t,u,interp)
  T = eltype(eltype(u))
  N = length((size(u[1])..., length(u)))
  TestSolution{T,N,true,typeof(t),typeof(u),typeof(interp)}(t,u,interp,true,:Success)
end
TestSolution(interp::DESolution) = TestSolution{Void,0,true,Void,Void,typeof(interp)}(nothing,nothing,interp,true,:Success)
hasinterp{T,N,hi,tType,uType,iType}(::TestSolution{T,N,hi,tType,uType,iType}) = hi
"""
`appxtrue(sol::AbstractODESolution,sol2::TestSolution)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.
"""
function appxtrue(sol::AbstractODESolution,sol2::TestSolution)
  if sol2.u == nothing && hasinterp(sol2)
    _sol = TestSolution(sol.t,sol2(sol.t),sol2)
  else
    _sol = sol2
  end

  errors = Dict(:final=>recursive_mean(abs.(sol[end]-_sol[end])))
  if _sol.dense
    timeseries_analytic = _sol(sol.t)
    errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol-timeseries_analytic))
    errors[:l2] = sqrt(recursive_mean(vecvecapply((x)->float(x).^2,sol-timeseries_analytic)))
    densetimes = collect(linspace(sol.t[1],sol.t[end],100))
    interp_u = sol(densetimes)
    interp_analytic = _sol(densetimes)
    interp_errors = Dict(:L∞=>maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic)),
                         :L2=>sqrt(recursive_mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic))))
    errors = merge(errors,interp_errors)
  else
    timeseries_analytic = sol2.u
    if sol.t == sol2.t
      errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol-timeseries_analytic))
      errors[:l2] = sqrt(recursive_mean(vecvecapply((x)->float(x).^2,sol-timeseries_analytic)))
    end
  end
  build_solution(sol,timeseries_analytic,errors)
end

function appxtrue(sol::AbstractFEMSolution,sol2::AbstractFEMSolution)
  u_analytic = sol2[end]
  errors = Dict(:l∞=>maximum(abs.(sol[end]-u_analytic)),:l2=>norm(sol[end]-u_analytic,2))
  FEMSolution(sol,u_analytic,errors)
end

"""
`appxtrue(sol::AbstractODESolution,sol2::AbstractODESolution)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.
"""
function appxtrue(sol::AbstractODESolution,sol2::AbstractODESolution;timeseries_errors=sol2.dense,dense_errors=sol2.dense)
  errors = Dict(:final=>recursive_mean(abs.(sol[end]-sol2[end])))
  if sol2.dense
    timeseries_analytic = sol2(sol.t)
    errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol-timeseries_analytic))
    errors[:l2] = sqrt(recursive_mean(vecvecapply((x)->float(x).^2,sol-timeseries_analytic)))
    if dense_errors
      densetimes = collect(linspace(sol.t[1],sol.t[end],100))
      interp_u = sol(densetimes)
      interp_analytic = sol2(densetimes)
      interp_errors = Dict(:L∞=>maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic)),
                           :L2=>sqrt(recursive_mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic))))
      errors = merge(errors,interp_errors)
    end
  else
    timeseries_analytic = sol2.u
    if timeseries_errors && sol.t == sol2.t
      errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol-timeseries_analytic))
      errors[:l2] = sqrt(recursive_mean(vecvecapply((x)->float(x).^2,sol-timeseries_analytic)))
    end
  end
  build_solution(sol,timeseries_analytic,errors)
end

function appxtrue(sim::MonteCarloSolution,appx_setup;kwargs...)
  _new_sols = Vector{DESolution}(length(sim.u))
  for i in eachindex(sim)
    @show i
    prob = sim[i].prob
    prob2 = SDEProblem(prob.f,prob.g,prob.u0,prob.tspan,noise=NoiseWrapper(sim[i].W))
    true_sol = solve(prob2,appx_setup[:alg];appx_setup...)
    _new_sols[i] = appxtrue(sim[i],true_sol)
  end
  new_sols = convert(Vector{typeof(_new_sols[1])},_new_sols)
  calculate_monte_errors(new_sols;converged=sim.converged,elapsedTime=sim.elapsedTime,kwargs...)
end
