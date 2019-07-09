using Statistics
## Shootouts

mutable struct Shootout
  setups::Vector{Dict{Symbol,Any}}
  times#::Vector{Float64}
  errors#::Vector{uType}
  effs#::Vector{Float64} # Efficiencies
  effratios#::Matrix{uEltype}
  solutions
  names::Vector{String}
  N::Int
  bestidx::Int
  winner::String
end

mutable struct ShootoutSet
  shootouts::Vector{Shootout}
  probs#::Vector{DEProblem}
  probaux#::Vector{Dict{Symbol,Any}}
  N::Int
  winners::Vector{String}
end

function ode_shootout(args...;kwargs...)
  @warn("ode_shootout is deprecated. Use ShootOut instead")
  ShootOut(args...;kwargs...)
end

function Shootout(prob,setups;appxsol=nothing,names=nothing,error_estimate=:final,numruns=20,seconds=2,kwargs...)
  N = length(setups)
  errors = Vector{Float64}(undef,N)
  solutions = Vector{Any}(undef,N)
  effs = Vector{Float64}(undef,N)
  times = Vector{Float64}(undef,N)
  effratios = Matrix{Float64}(undef,N,N)
  timeseries_errors = error_estimate ∈ TIMESERIES_ERRORS
  dense_errors = error_estimate ∈ DENSE_ERRORS
  if names == nothing
    names = [string(nameof(typeof(setup[:alg]))) for setup in setups]
  end
  for i in eachindex(setups)
    sol = solve(prob,setups[i][:alg];timeseries_errors=timeseries_errors,
    dense_errors = dense_errors,kwargs...,setups[i]...) # Compile and get result

    if :prob_choice ∈ keys(setups[i])
      cur_appxsol = appxsol[setups[i][:prob_choice]]
    else
      cur_appxsol = appxsol
    end

    if cur_appxsol != cur_appxsol
      errsol = appxtrue(sol,cur_appxsol)
      errors[i] = errsol.errors[error_estimate]
      solutions[i] = errsol
    else
      errors[i] = sol.errors[error_estimate]
      solutions[i] = sol
    end

    if haskey(setups[i], :prob_choice)
      _prob = prob[setups[i][:prob_choice]]
    else
      _prob = prob
    end

    benchmark_f = let _prob=_prob,alg=setups[i][:alg],sol=sol,kwargs=kwargs
      function benchmark_f()
        @elapsed solve(_prob,alg,(sol.u),(sol.t),(sol.k);
              timeseries_errors = false,
              dense_errors = false, kwargs...)
      end
    end

    b_t =  benchmark_f()
    if b_t > seconds
      times[i] = b_t
    else
      times[i] = minimum([b_t;map(i->benchmark_f(),2:numruns)])
    end

    effs[i] = 1/(errors[i]*times[i])
  end
  for j in 1:N, i in 1:N
    effratios[i,j] = effs[i]/effs[j]
  end
  bestidx = findall((y)->y==maximum(effs),effs)[1]; winner = names[bestidx]
  return Shootout(setups,times,errors,effs,effratios,solutions,names,N,bestidx,winner)
end

function ode_shootoutset(args...;kwargs...)
  @warn("ode_shootoutset is deprecated. Use ShootoutSet instead")
  ShootoutSet(args...;kwargs...)
end

function ShootoutSet(probs,setups;probaux=nothing,
                     names=nothing,print_names=false,kwargs...)
  N = length(probs)
  shootouts = Vector{Shootout}(undef,N)
  winners = Vector{String}(undef,N)
  if names == nothing
    names = [string(nameof(typeof(setup[:alg]))) for setup in setups]
  end
  if probaux == nothing
    probaux = Vector{Dict{Symbol,Any}}(undef,N)
    for i in 1:N
      probaux[i] = Dict{Symbol,Any}()
    end
  end
  for i in eachindex(probs)
    print_names && println(names[i])
    shootouts[i] = Shootout(probs[i],setups;names=names,kwargs...,probaux[i]...)
    winners[i] = shootouts[i].winner
  end
  return ShootoutSet(shootouts,probs,probaux,N,winners)
end

Base.length(shoot::Shootout) = shoot.N
Base.size(shoot::Shootout) = length(shoot)
Base.getindex(shoot::Shootout,i::Int) = shoot.effs[i]
Base.getindex(shoot::Shootout,::Colon) = shoot.effs
Base.firstindex(shoot::Shootout) = 1
Base.lastindex(shoot::Shootout) = lastindex(shoot.effs)

function Base.show(io::IO, shoot::Shootout)
  println(io,"Winner: $(shoot.winner)")
  println(io,"EffRatios: $(shoot.effratios[shoot.bestidx,:])")
end

Base.length(set::ShootoutSet) = set.N
Base.size(set::ShootoutSet) = length(set)
Base.getindex(set::ShootoutSet,i::Int) = set.shootouts[i]
Base.getindex(set::ShootoutSet,::Colon) = set.shootouts
Base.show(io::IO, set::ShootoutSet) = print(io,"ShootoutSet of $(set.N) shootouts ")
Base.firstindex(shoot::ShootoutSet) = 1
Base.lastindex(shoot::ShootoutSet) = lastindex(shoot.shootouts)

## WorkPrecisions

mutable struct WorkPrecision
  prob
  abstols
  reltols
  errors
  times
  name
  N::Int
end

mutable struct WorkPrecisionSet
  wps::Vector{WorkPrecision}
  N::Int
  abstols
  reltols
  prob
  setups
  names
  sample_error
  error_estimate
  numruns
end

function WorkPrecision(prob,alg,abstols,reltols,dts=nothing;
                       name=nothing,appxsol=nothing,error_estimate=:final,numruns=20,seconds=2,kwargs...)
  N = length(abstols)
  errors = Vector{Float64}(undef,N)
  times = Vector{Float64}(undef,N)
  if name == nothing
    name = "WP-Alg"
  end

  if :prob_choice ∈ keys(kwargs)
    _prob = prob[kwargs[:prob_choice]]
  else
    _prob = prob
  end

  let _prob = _prob
    timeseries_errors = error_estimate ∈ TIMESERIES_ERRORS
    dense_errors = error_estimate ∈ DENSE_ERRORS
    for i in 1:N
      # Calculate errors and precompile
      if dts == nothing
        sol = solve(_prob,alg;kwargs...,abstol=abstols[i],
        reltol=reltols[i],timeseries_errors=timeseries_errors,
        dense_errors = dense_errors) # Compile and get result
      else
        sol = solve(_prob,alg;kwargs...,abstol=abstols[i],
        reltol=reltols[i],dt=dts[i],timeseries_errors=timeseries_errors,
        dense_errors = dense_errors) # Compile and get result
      end

      if haskey(kwargs, :prob_choice)
        cur_appxsol = appxsol[kwargs[:prob_choice]]
      else
        cur_appxsol = appxsol
      end

      if cur_appxsol != nothing
        errsol = appxtrue(sol,cur_appxsol)
        errors[i] = mean(errsol.errors[error_estimate])
      else
        errors[i] = mean(sol.errors[error_estimate])
      end

      benchmark_f = let dts=dts,_prob=_prob,alg=alg,sol=sol,abstols=abstols,reltols=reltols,kwargs=kwargs
        function benchmark_f()
          if dts == nothing
            @elapsed solve(_prob,alg,(sol.u),(sol.t),(sol.k);
                  abstol=(abstols[i]),
                  reltol=(reltols[i]),
                  timeseries_errors = false,
                  dense_errors = false, kwargs...)
          else
            @elapsed solve(_prob,alg,(sol.u),(sol.t),(sol.k);
                  abstol=(abstols[i]),
                  reltol=(reltols[i]),
                  dt=(dts[i]),
                  timeseries_errors = false,
                  dense_errors = false, kwargs...)
          end
        end
      end

      b_t =  benchmark_f()
      if b_t > seconds
        times[i] = b_t
      else
        times[i] = minimum([b_t;map(i->benchmark_f(),2:numruns)])
      end
    end
  end
  return WorkPrecision(prob,abstols,reltols,errors,times,name,N)
end

function WorkPrecisionSet(prob,
                          abstols,reltols,setups;
                          print_names=false,names=nothing,appxsol=nothing,
                          error_estimate=:final,
                          test_dt=nothing,kwargs...)
  N = length(setups)
  wps = Vector{WorkPrecision}(undef,N)
  if names == nothing
    names = [string(nameof(typeof(setup[:alg]))) for setup in setups]
  end
  for i in 1:N
    print_names && println(names[i])
    if haskey(setups[i],:dts)
      wps[i] = WorkPrecision(prob,setups[i][:alg],abstols,reltols,setups[i][:dts];
                                 appxsol=appxsol,
                                 error_estimate=error_estimate,
                                 name=names[i],kwargs...,setups[i]...)
    else
      wps[i] = WorkPrecision(prob,setups[i][:alg],abstols,reltols;
                                 appxsol=appxsol,
                                 error_estimate=error_estimate,
                                 name=names[i],kwargs...,setups[i]...)
    end
  end
  return WorkPrecisionSet(wps,N,abstols,reltols,prob,setups,names,nothing,error_estimate,nothing)
end

@def error_calculation begin
  if !DiffEqBase.has_analytic(prob.f)
    t = prob.tspan[1]:test_dt:prob.tspan[2]
    brownian_values = cumsum([[zeros(size(prob.u0))];[sqrt(test_dt)*randn(size(prob.u0)) for i in 1:length(t)-1]])
    brownian_values2 = cumsum([[zeros(size(prob.u0))];[sqrt(test_dt)*randn(size(prob.u0)) for i in 1:length(t)-1]])
    np = NoiseGrid(t,brownian_values,brownian_values2)
    _prob = remake(prob,noise=np);
    true_sol = solve(_prob,appxsol_setup[:alg];kwargs...,appxsol_setup...)
  else
    _prob = prob
  end

  # Get a cache
  if !haskey(setups[1],:dts)
    sol = solve(_prob,setups[1][:alg];
          kwargs...,setups[1]...,
          abstol=abstols[1],
          reltol=reltols[1],
          timeseries_errors=false,
          dense_errors = false)
  else
    sol = solve(_prob,setups[1][:alg];
          kwargs...,setups[1]...,abstol=abstols[1],
          reltol=reltols[1],dt=setups[1][:dts][1],
          timeseries_errors=false,
          dense_errors = false)
  end

  for j in 1:M, k in 1:N
    if !haskey(setups[k],:dts)
      sol = solve(_prob,setups[k][:alg];
            kwargs...,setups[k]...,
            abstol=abstols[j],
            reltol=reltols[j],
            timeseries_errors=timeseries_errors,
            dense_errors = dense_errors)
    else
      sol = solve(_prob,setups[k][:alg];
            kwargs...,setups[k]...,abstol=abstols[j],
            reltol=reltols[j],dt=setups[k][:dts][j],
            timeseries_errors=timeseries_errors,
            dense_errors = dense_errors)
    end
    DiffEqBase.has_analytic(prob.f) ? err_sol = sol : err_sol = appxtrue(sol,true_sol)
    tmp_solutions[i,j,k] = err_sol
  end
end

function WorkPrecisionSet(prob::AbstractRODEProblem,abstols,reltols,setups,test_dt=nothing;
                          numruns=20,numruns_error = 20,
                          print_names=false,names=nothing,appxsol_setup=nothing,
                          error_estimate=:final,parallel_type = :none,kwargs...)

  timeseries_errors = DiffEqBase.has_analytic(prob.f) && error_estimate ∈ TIMESERIES_ERRORS
  weak_timeseries_errors = error_estimate ∈ WEAK_TIMESERIES_ERRORS
  weak_dense_errors = error_estimate ∈ WEAK_DENSE_ERRORS
  dense_errors = DiffEqBase.has_analytic(prob.f) && error_estimate ∈ DENSE_ERRORS
  N = length(setups); M = length(abstols)
  times = Array{Float64}(undef,M,N)
  tmp_solutions = Array{Any}(undef,numruns_error,M,N)
  if names == nothing
    names = [string(nameof(typeof(setup[:alg]))) for setup in setups]
  end
  time_tmp = Vector{Float64}(undef,numruns)

  # First calculate all of the errors
  if parallel_type == :threads
    Threads.@threads for i in 1:numruns_error
      @error_calculation
    end
  elseif parallel_type == :none
    for i in 1:numruns_error
      @info "Error calculation: $i/$numruns_error"
      @error_calculation
    end
  end
  analytical_solution_ends = [tmp_solutions[i,1,1].u_analytic[end] for i in 1:numruns_error]
  sample_error = 1.96std(norm.(analytical_solution_ends))/sqrt(numruns_error)
  _solutions_k = [[EnsembleSolution(tmp_solutions[:,j,k],0.0,true) for j in 1:M] for k in 1:N]
  solutions = [[DiffEqBase.calculate_ensemble_errors(sim;weak_timeseries_errors=weak_timeseries_errors,weak_dense_errors=weak_dense_errors) for sim in sol_k] for sol_k in _solutions_k]
  if error_estimate ∈ WEAK_ERRORS
    errors = [[solutions[j][i].weak_errors[error_estimate] for i in 1:M] for j in 1:N]
  else
    errors = [[solutions[j][i].error_means[error_estimate] for i in 1:M] for j in 1:N]
  end

  local _sol

  # Now time it
  for k in 1:N
    # precompile
    GC.gc()
    if !haskey(setups[k],:dts)
      _sol = solve(prob,setups[k][:alg];
            kwargs...,setups[k]...,
            abstol=abstols[1],
            reltol=reltols[1],
            timeseries_errors=false,
            dense_errors = false)
    else
      _sol = solve(prob,setups[k][:alg];
            kwargs...,setups[k]...,abstol=abstols[1],
            reltol=reltols[1],dt=setups[k][:dts][1],
            timeseries_errors=false,
            dense_errors = false)
    end
    x = isempty(_sol.t) ? 0 : round(Int,mean(_sol.t) - sum(_sol.t)/length(_sol.t))
    GC.gc()
    for j in 1:M
      for i in 1:numruns
        time_tmp[i] = @elapsed if !haskey(setups[k],:dts)
          sol = solve(prob,setups[k][:alg];
                kwargs...,setups[k]...,
                abstol=abstols[j],
                reltol=reltols[j],
                timeseries_errors=false,
                dense_errors = false)
        else
          sol = solve(prob,setups[k][:alg];
                kwargs...,setups[k]...,abstol=abstols[j],
                reltol=reltols[j],dt=setups[k][:dts][j],
                timeseries_errors=false,
                dense_errors = false)
        end
      end
      times[j,k] = mean(time_tmp) + x
      GC.gc()
    end
  end

  wps = [WorkPrecision(prob,abstols,reltols,errors[i],times[:,i],names[i],N) for i in 1:N]
  WorkPrecisionSet(wps,N,abstols,reltols,prob,setups,names,sample_error,error_estimate,numruns_error)
end

@def sample_errors begin
  if !DiffEqBase.has_analytic(prob.f)
    true_sol = solve(prob,appxsol_setup[:alg];kwargs...,appxsol_setup...,
                     save_everystep=false)
    analytical_solution_ends[i] = norm(true_sol.u[end])
  else
    _dt = prob.tspan[2] - prob.tspan[1]
    if typeof(prob.u0) <: Number
      W = sqrt(_dt)*randn()
    else
      W = sqrt(_dt)*randn(size(prob.u0))
    end
    analytical_solution_ends[i] = norm(prob.f.analytic(prob.u0,prob.p,prob.tspan[2],W))
  end
end

function get_sample_errors(prob::AbstractRODEProblem,test_dt=nothing;
                          appxsol_setup=nothing,
                          numruns=20,std_estimation_runs = maximum(numruns),
                          error_estimate=:final,parallel_type = :none,kwargs...)
  _std_estimation_runs = Int(std_estimation_runs)
  analytical_solution_ends = Vector{typeof(norm(prob.u0))}(undef,_std_estimation_runs)
  if parallel_type == :threads
    Threads.@threads for i in 1:_std_estimation_runs
      @sample_errors
    end
  elseif parallel_type == :none
    for i in 1:_std_estimation_runs
      @info "Standard deviation estimation: $i/$_std_estimation_runs"
      @sample_errors
    end
  end
  est_std = std(analytical_solution_ends)
  if typeof(numruns) <: Number
    return 1.96est_std/sqrt(numruns)
  else
    return [1.96est_std/sqrt(_numruns) for _numruns in numruns]
  end
end

Base.length(wp::WorkPrecision) = wp.N
Base.size(wp::WorkPrecision) = length(wp)
Base.getindex(wp::WorkPrecision,i::Int) = wp.times[i]
Base.getindex(wp::WorkPrecision,::Colon) = wp.times
Base.firstindex(wp::WorkPrecision) = 1
Base.lastindex(wp::WorkPrecision) = lastindex(wp.times)

function Base.show(io::IO, wp::WorkPrecision)
  println(io,"Name: $(wp.name)")
  println(io,"Times: $(wp.times)")
  println(io,"Errors: $(wp.errors)")
end

Base.length(wp_set::WorkPrecisionSet) = wp_set.N
Base.size(wp_set::WorkPrecisionSet) = length(wp_set)
Base.getindex(wp_set::WorkPrecisionSet,i::Int) = wp_set.wps[i]
Base.getindex(wp_set::WorkPrecisionSet,::Colon) = wp_set.wps
Base.firstindex(wp_set::WorkPrecisionSet) = 1
Base.lastindex(wp_set::WorkPrecisionSet)  = lastindex(wp_set.wps)

Base.show(io::IO, wp_set::WorkPrecisionSet) = println(io,"WorkPrecisionSet of $(wp_set.N) wps")
