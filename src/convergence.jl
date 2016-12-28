type ConvergenceSimulation{SolType<:DESolution}
  solutions::Array{SolType}
  errors
  N
  auxdata
  𝒪est
  convergence_axis
end

function ConvergenceSimulation(solutions,convergence_axis;auxdata=nothing,additional_errors=nothing)
  N = size(solutions,1)
  uEltype = eltype(solutions[1].u[1])
  errors = Dict() #Should add type information
  for k in keys(solutions[1].errors)
    errors[k] = reshape(uEltype[sol.errors[k] for sol in solutions],size(solutions)...)
  end
  if additional_errors != nothing
    for k in keys(additional_errors)
      errors[k] = additional_errors[k]
    end
  end
  𝒪est = Dict(map(calc𝒪estimates,errors))
  𝒪esttmp = Dict() #Makes Dict of Any to be more compatible
  for (k,v) in 𝒪est
    if length(v)==1 push!(𝒪esttmp,Pair(k,v[1]))
    else push!(𝒪esttmp,Pair(k,v))
    end
  end
  𝒪est = 𝒪esttmp
  return(ConvergenceSimulation(solutions,errors,N,auxdata,𝒪est,convergence_axis))
end

function test_convergence(dts::AbstractArray,prob::AbstractSDEProblem,alg;numMonte=10000,save_timeseries=true,timeseries_steps=1,timeseries_errors=save_timeseries,adaptive=false,kwargs...)
  N = length(dts)
  is = repmat(1:N,1,numMonte)'
  solutions = pmap((i)->solve(prob,alg;dt=dts[i],save_timeseries=save_timeseries,timeseries_steps=timeseries_steps,adaptive=adaptive,timeseries_errors=timeseries_errors,kwargs...),is)
  if typeof(prob) <: SDEProblem
    solutions = convert(Array{SDESolution},solutions)
  elseif typeof(prob) <: SDETestProblem
    solutions = convert(Array{SDETestSolution},solutions)
  end
  solutions = reshape(solutions,numMonte,N)
  auxdata = Dict("dts" =>  dts)
  # Now Calculate Weak Errors
  additional_errors = Dict()
  # Final
  m_final = mean([s[end] for s in solutions],1)
  m_final_analytic = mean([s.u_analytic[end] for s in solutions],1)
  additional_errors[:weak_final] = mean.(vecvecapply((x)->abs.(x),m_final - m_final_analytic))
  if timeseries_errors
    l2_tmp = Vector{eltype(solutions[1][1])}(size(solutions,2))
    max_tmp = Vector{eltype(solutions[1][1])}(size(solutions,2))
    for i in 1:size(solutions,2)
      solcol = @view solutions[:,i]
      m_errors = [mean([solcol[j][i] for j in 1:length(solcol)]) for i in 1:length(solcol[1])]
      m_errors_analytic = [mean([solcol[j].u_analytic[i] for j in 1:length(solcol)]) for i in 1:length(solcol[1])]
      ts_weak_errors = [abs.(m_errors[i] - m_errors_analytic[i]) for i in 1:length(m_errors)]
      ts_l2_errors = [sqrt.(sum(abs2,err)/length(err)) for err in ts_weak_errors]
      l2_tmp[i] = sqrt(sum(abs2,ts_l2_errors)/length(ts_l2_errors))
      max_tmp[i] = maximum([maximum(err) for err in ts_weak_errors])
    end
    additional_errors[:weak_l2] = l2_tmp
    additional_errors[:weak_l∞] = max_tmp
  end
  ConvergenceSimulation(solutions,dts,auxdata=auxdata,additional_errors=additional_errors)
end

function test_convergence(dts::AbstractArray,prob::AbstractODEProblem,alg;save_timeseries=true,adaptive=false,kwargs...)
  N = length(dts)
  solutions = [solve(prob,alg;dt=dts[i],save_timeseries=save_timeseries,adaptive=adaptive,kwargs...) for i=1:N]
  auxdata = Dict(:dts =>  dts)
  ConvergenceSimulation(solutions,dts,auxdata=auxdata)
end

#=
function test_convergence(dts::AbstractArray,dxs::AbstractArray,prob::AbstractHeatProblem,convergence_axis;T=1,alg=:Euler)
  if length(dts)!=length(dxs) error("Lengths of dts!=dxs. Invalid convergence simulation") end
  solutions = [solve(parabolic_squaremesh([0 1 0 1],dxs[i],dts[i],T,:dirichlet),prob,alg=alg) for i in eachindex(dts)]
  auxdata = Dict(
            :dts => [sol.fem_mesh.dt for sol in solutions],
            :dxs => [sol.fem_mesh.dx for sol in solutions],
            :Δμs => [sol.fem_mesh.μ  for sol in solutions],
            :Δνs => [sol.fem_mesh.ν  for sol in solutions])
  return(ConvergenceSimulation(solutions,convergence_axis,auxdata=auxdata))
end
=#

function test_convergence(probs,convergence_axis;kwargs...)
  ConvergenceSimulation([solve(prob;kwargs...) for prob in probs],convergence_axis)
end

function test_convergence(c::ConvergenceSetup;kwargs...)
  test_convergence(c.probs,c.convergence_axis;kwargs...)
end

function calc𝒪estimates(error::Pair)
  key = error.first
  error =error.second
  if ndims(error)>1 error=mean(error,1) end
  S = Vector{eltype(error)}(length(error)-1)
  for i=1:length(error)-1
    S[i] = log2(error[i+1]/error[i])
  end
  return(Pair(key,abs.(mean(S,1))))
end

"""
length(simres::ConvergenceSimulation)

Returns the number of simultations in the Convergence Simulation
"""
Base.length(sim::ConvergenceSimulation) = sim.N
Base.endof( sim::ConvergenceSimulation) = length(sim)
Base.getindex(sim::ConvergenceSimulation,i::Int) = sim.solutions[i]
Base.getindex(sim::ConvergenceSimulation,i::Int,I::Int...) = sim.solutions[i][I]

Base.length(sim::ConvergenceSetup) = sim.probs
