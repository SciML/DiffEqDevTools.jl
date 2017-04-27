"""
TestSolution

"""
type TestSolution{T,N,hasinterp,tType,uType,iType} <: AbstractTimeseriesSolution{T,N}
  t::tType
  u::uType
  interp::iType
  dense::Bool
end
(T::TestSolution)(t) = T.interp(t)
function TestSolution(t,u)
  T = eltype(eltype(u))
  N = length((size(u[1])..., length(u)))
  TestSolution{T,N,false,typeof(t),typeof(u),Void}(t,u,nothing,false)
end
function TestSolution(t,u,interp)
  T = eltype(eltype(u))
  N = length((size(u[1])..., length(u)))
  TestSolution{T,N,length(true),typeof(t),typeof(u),typeof(interp)}(t,u,interp,true)
end
TestSolution(interp::DESolution) = TestSolution{Void,0,true,Void,Void,typeof(interp)}(nothing,nothing,interp,true)

"""
`appxtrue(sol::AbstractODESolution,sol2::TestSolution)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.
"""
function appxtrue{hasinterp}(sol::AbstractODESolution,sol2::TestSolution{hasinterp})
  if sol2.u == nothing && hasinterp
    sol2.u = sol2(sol.t)
    sol2.t = sol.t
  end
  errors = Dict(:final=>recursive_mean(abs.(sol[end]-sol2[end])))
  if sol2.dense
    timeseries_analytic = sol2(sol.t)
    errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol[:]-timeseries_analytic))
    errors[:l2] = sqrt(recursive_mean(vecvecapply((x)->float(x).^2,sol[:]-timeseries_analytic)))
    if !(typeof(sol) <: AbstractRODESolution) && sol.dense
      densetimes = collect(linspace(sol.t[1],sol.t[end],100))
      interp_u = sol(densetimes)
      interp_analytic = sol2(densetimes)
      interp_errors = Dict(:L∞=>maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic)),
                           :L2=>sqrt(recursive_mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic))))
      errors = merge(errors,interp_errors)
    end
  end
  build_solution(sol,sol2.u,errors)
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
function appxtrue(sol::AbstractODESolution,sol2::AbstractODESolution)
  errors = Dict(:final=>recursive_mean(abs.(sol[end]-sol2[end])))
  if !(typeof(sol2) <: AbstractRODESolution) && sol2.dense
    timeseries_analytic = sol2(sol.t)
    errors = Dict(:final=>recursive_mean(abs.(sol[end]-sol2[end])),:l∞=>maximum(vecvecapply((x)->abs.(x),sol[:]-timeseries_analytic)),:l2=>sqrt(recursive_mean(vecvecapply((x)->float(x).^2,sol[:]-timeseries_analytic))))
    if !(typeof(sol) <: AbstractRODESolution) && sol.dense
      densetimes = collect(linspace(sol.t[1],sol.t[end],100))
      interp_u = sol(densetimes)
      interp_analytic = sol2(densetimes)
      interp_errors = Dict(:L∞=>maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic)),:L2=>sqrt(recursive_mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic))))
      errors = merge(errors,interp_errors)
    end
  end
  build_solution(sol,sol2.u,errors)
end
