type ConvergenceSimulation{SolType}
  solutions::Array{SolType}
  errors
  N
  auxdata
  𝒪est
  convergence_axis
end

function ConvergenceSimulation(solutions,convergence_axis;
                               auxdata=nothing,additional_errors=nothing)
  N = size(solutions,1)
  uEltype = eltype(solutions[1].u[1])
  errors = Dict() #Should add type information
  if isempty(solutions[1].errors)
    error("Errors dictionary is empty. No analytical solution set.")
  end
  for k in keys(solutions[1].errors)
    errors[k] = [mean(sol.errors[k]) for sol in solutions]
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

function test_convergence(dts::AbstractArray,prob::Union{AbstractRODEProblem,AbstractSDEProblem},
                          alg;numMonte=10000,save_everystep=true,timeseries_steps=1,
                          timeseries_errors=save_everystep,adaptive=false,
                          weak_timeseries_errors=false,weak_dense_errors=false,kwargs...)
  N = length(dts)
  monte_prob = MonteCarloProblem(prob)
  _solutions = [solve(monte_prob,alg;dt=dts[i],save_everystep=save_everystep,
                      timeseries_steps=timeseries_steps,adaptive=adaptive,
                      timeseries_errors=timeseries_errors,num_monte=numMonte,
                      kwargs...) for i in 1:N]
  solutions = [calculate_monte_errors(sim;weak_timeseries_errors=weak_timeseries_errors,weak_dense_errors=weak_dense_errors) for sim in _solutions]
  auxdata = Dict("dts" =>  dts)
  # Now Calculate Weak Errors
  additional_errors = Dict()
  for k in keys(solutions[1].weak_errors)
    additional_errors[k] = [sol.weak_errors[k] for sol in solutions]
  end
  ConvergenceSimulation(solutions,dts,auxdata=auxdata,additional_errors=additional_errors)
end

function test_convergence(dts::AbstractArray,prob::AbstractODEProblem,alg;save_everystep=true,adaptive=false,kwargs...)
  N = length(dts)
  solutions = [solve(prob,alg;dt=dts[i],save_everystep=save_everystep,adaptive=adaptive,kwargs...) for i=1:N]
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

function test_convergence(probs,convergence_axis,alg;kwargs...)
  ConvergenceSimulation([solve(prob,alg;kwargs...) for prob in probs],convergence_axis)
end

function test_convergence(c::ConvergenceSetup,alg::DEAlgorithm;kwargs...)
  test_convergence(c.probs,c.convergence_axis,alg;kwargs...)
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
