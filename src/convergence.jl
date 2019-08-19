mutable struct ConvergenceSimulation{SolType}
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
  𝒪est = Dict((calc𝒪estimates(p) for p = pairs(errors)))
  #𝒪est = Dict(map(calc𝒪estimates,errors))
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
                          alg;trajectories,save_everystep=true,timeseries_steps=1,
                          timeseries_errors=save_everystep,adaptive=false,
                          weak_timeseries_errors=false,weak_dense_errors=false,kwargs...)
  N = length(dts)
  ensemble_prob = EnsembleProblem(prob)
  _solutions = [solve(ensemble_prob,alg;dt=dts[i],save_everystep=save_everystep,
                      timeseries_steps=timeseries_steps,adaptive=adaptive,
                      timeseries_errors=timeseries_errors,trajectories=trajectories,
                      kwargs...) for i in 1:N]
  solutions = [DiffEqBase.calculate_ensemble_errors(sim;weak_timeseries_errors=weak_timeseries_errors,weak_dense_errors=weak_dense_errors) for sim in _solutions]
  auxdata = Dict("dts" =>  dts)
  # Now Calculate Weak Errors
  additional_errors = Dict()
  for k in keys(solutions[1].weak_errors)
    additional_errors[k] = [sol.weak_errors[k] for sol in solutions]
  end
  ConvergenceSimulation(solutions,dts,auxdata=auxdata,additional_errors=additional_errors)
end

function analyticless_test_convergence(dts::AbstractArray,
                          prob::Union{AbstractRODEProblem,AbstractSDEProblem},
                          alg,test_dt;trajectories=100,
                          save_everystep=true,timeseries_steps=1,
                          timeseries_errors=save_everystep,adaptive=false,
                          weak_timeseries_errors=false,weak_dense_errors=false,kwargs...)
  _solutions = []
  tmp_solutions = Array{Any}(undef,trajectories,length(dts))
  for j in 1:trajectories
    @info "Monte Carlo iteration: $j/$trajectories"
    t = prob.tspan[1]:test_dt:prob.tspan[2]
    if prob.noise_rate_prototype === nothing
      brownian_values = cumsum([[zeros(size(prob.u0))];[sqrt(test_dt)*randn(size(prob.u0)) for i in 1:length(t)-1]])
      brownian_values2 = cumsum([[zeros(size(prob.u0))];[sqrt(test_dt)*randn(size(prob.u0)) for i in 1:length(t)-1]])
    else
      brownian_values = cumsum([[zeros(size(prob.noise_rate_prototype,2))];[sqrt(test_dt)*randn(size(prob.noise_rate_prototype,2)) for i in 1:length(t)-1]])
      brownian_values2 = cumsum([[zeros(size(prob.noise_rate_prototype,2))];[sqrt(test_dt)*randn(size(prob.noise_rate_prototype,2)) for i in 1:length(t)-1]])
    end
    np = NoiseGrid(t,brownian_values,brownian_values2)
    _prob = SDEProblem(prob.f,prob.g,prob.u0,prob.tspan,
                       noise=np,
                       noise_rate_prototype=prob.noise_rate_prototype);
    true_sol = solve(_prob,alg;adaptive=adaptive,dt=test_dt);
    for i in 1:length(dts)
      sol = solve(_prob,alg;dt=dts[i],adaptive=adaptive);
      err_sol = appxtrue(sol,true_sol)
      tmp_solutions[j,i] = err_sol
    end
  end
  _solutions = [EnsembleSolution(tmp_solutions[:,i],0.0,true) for i in 1:length(dts)]
  solutions = [DiffEqBase.calculate_ensemble_errors(sim;weak_timeseries_errors=weak_timeseries_errors,weak_dense_errors=weak_dense_errors) for sim in _solutions]
  auxdata = Dict("dts" =>  dts)
  # Now Calculate Weak Errors
  additional_errors = Dict()
  for k in keys(solutions[1].weak_errors)
    additional_errors[k] = [sol.weak_errors[k] for sol in solutions]
  end
  ConvergenceSimulation(solutions,dts,auxdata=auxdata,additional_errors=additional_errors)
end

function test_convergence(dts::AbstractArray,prob::AbstractODEProblem,alg;
                          save_everystep=true,adaptive=false,kwargs...)
  N = length(dts)
  solutions = [solve(prob,alg;dt=dts[i],save_everystep=save_everystep,adaptive=adaptive,kwargs...) for i=1:N]
  auxdata = Dict(:dts =>  dts)
  ConvergenceSimulation(solutions,dts,auxdata=auxdata)
end

function analyticless_test_convergence(dts::AbstractArray,prob::AbstractODEProblem,
                                       alg,appxsol_setup;
                                       save_everystep=true,adaptive=false,kwargs...)
  true_sol = solve(prob,appxsol_setup[:alg];appxsol_setup...);
  N = length(dts)
  _solutions = [solve(prob,alg;dt=dts[i],save_everystep=save_everystep,adaptive=adaptive,kwargs...) for i=1:N]
  solutions = [appxtrue(sol,true_sol) for sol in _solutions]
  auxdata = Dict(:dts =>  dts)
  ConvergenceSimulation(solutions,dts,auxdata=auxdata)
end

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
  S = Vector{eltype(error)}(undef, length(error)-1)
  for i=1:length(error)-1
    S[i] = log2(error[i+1]/error[i])
  end
  return(Pair(key,abs.(mean(S,dims=1))))
end

"""
length(simres::ConvergenceSimulation)

Returns the number of simultations in the Convergence Simulation
"""
Base.length(sim::ConvergenceSimulation) = sim.N
Base.getindex(sim::ConvergenceSimulation,i::Int) = sim.solutions[i]
Base.getindex(sim::ConvergenceSimulation,i::Int,I::Int...) = sim.solutions[i][I]
Base.lastindex(sim::ConvergenceSimulation) = lastindex(sim.solutions)
Base.firstindex(sim::ConvergenceSimulation) = firstindex(sim.solutions)
Base.length(sim::ConvergenceSetup) = sim.probs
