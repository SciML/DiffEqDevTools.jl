type ConvergenceSimulation{SolType<:DESolution}
  solutions::Array{SolType}
  errors
  N
  auxdata
  𝒪est
  convergence_axis
end

function ConvergenceSimulation(solutions,convergence_axis;auxdata=nothing)
  N = size(solutions,1)
  uEltype = eltype(solutions[1].u[1])
  errors = Dict() #Should add type information
  for k in keys(solutions[1].errors)
    errors[k] = reshape(uEltype[sol.errors[k] for sol in solutions],size(solutions)...)
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

function test_convergence(dts::AbstractArray,prob::AbstractSDEProblem,alg;numMonte=10000,save_timeseries=true,timeseries_steps=1,adaptive=false,kwargs...)
  N = length(dts)
  is = repmat(1:N,1,numMonte)'
  solutions = pmap((i)->solve(prob,alg;dt=dts[i],save_timeseries=save_timeseries,timeseries_steps=timeseries_steps,adaptive=adaptive,kwargs...),is)
  if typeof(prob) <: SDEProblem
    solutions = convert(Array{SDESolution},solutions)
  elseif typeof(prob) <: SDETestProblem
    solutions = convert(Array{SDETestSolution},solutions)
  end
  solutions = reshape(solutions,numMonte,N)
  auxdata = Dict("dts" =>  dts)
  ConvergenceSimulation(solutions,dts,auxdata=auxdata)
end

function test_convergence(dts::AbstractArray,prob::AbstractODEProblem,alg;save_timeseries=true,adaptive=false,kwargs...)
  N = length(dts)
  solutions = [solve(prob,alg;dt=dts[i],save_timeseries=save_timeseries,adaptive=adaptive,kwargs...) for i=1:N]
  auxdata = Dict(:dts =>  dts)
  ConvergenceSimulation(solutions,dts,auxdata=auxdata)
end

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

function test_convergence(dxs::AbstractArray,prob::AbstractPoissonProblem)
  solutions = [solve(notime_squaremesh([0 1 0 1],dxs[i],:dirichlet),prob,solver=:Direct) for i in eachindex(dxs)]
  auxdata = Dict("dxs" => [sol.fem_mesh.dx for sol in solutions])
  return(ConvergenceSimulation(solutions,dxs,auxdata=auxdata))
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

function print(io::IO, sim::ConvergenceSimulation)
  println(io,"$(typeof(sim)) of length $(length(sim)).")
  print(io,"Convergence Estimates:")
  for (k,v) in sim.𝒪est
    print(" ($k,$v)")
  end
  println(io,"\n-----------Errors-----------")
  for (k,v) in sim.errors
    println(io,"$k: $v")
  end
end

function show(io::IO,sim::ConvergenceSimulation)
  println(io,"$(typeof(sim)) of length $(length(sim)).")
  print(io,"Convergence Estimates:")
  for (k,v) in sim.𝒪est
    print(io," ($k,$v)")
  end
end
