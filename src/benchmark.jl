## Shootouts

type Shootout
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

type ShootoutSet
  shootouts::Vector{Shootout}
  probs#::Vector{DEProblem}
  probaux#::Vector{Dict{Symbol,Any}}
  N::Int
  winners::Vector{String}
end

function ode_shootout(args...;kwargs...)
  warn("ode_shootout is deprecated. Use ShootOut instead")
  ShootOut(args...;kwargs...)
end

function Shootout(prob,setups;appxsol=nothing,numruns=20,names=nothing,error_estimate=:final,kwargs...)
  N = length(setups)
  errors = Vector{Float64}(N)
  solutions = Vector{DESolution}(N)
  effs = Vector{Float64}(N)
  times = Vector{Float64}(N)
  effratios = Matrix{Float64}(N,N)
  timeseries_errors = error_estimate ∈ TIMESERIES_ERRORS
  dense_errors = error_estimate ∈ DENSE_ERRORS
  if names == nothing
    names = [string(typeof(setups[i][:alg])) for i=1:N]
  end
  for i in eachindex(setups)
    sol = solve(prob,setups[i][:alg];timeseries_errors=timeseries_errors,
    dense_errors = dense_errors,kwargs...,setups[i]...) # Compile and get result
    sol = solve(prob,setups[i][:alg],sol.u,sol.t,sol.k;timeseries_errors=timeseries_errors,
    dense_errors = dense_errors,kwargs...,setups[i]...) # Compile and get result
    gc()
    t = @elapsed for j in 1:numruns
      sol = solve(prob,setups[i][:alg],sol.u,sol.t,sol.k;
                  kwargs...,setups[i]...,timeseries_errors=false,dense_errors=false)
    end
    if appxsol != nothing
      errsol = appxtrue(sol,appxsol)
      errors[i] = errsol.errors[error_estimate]
      solutions[i] = errsol
    else
      errors[i] = sol.errors[error_estimate]
      solutions[i] = sol
    end
    effs[i] = 1/(errors[i]*t)
    t = t/numruns
    times[i] = t
  end
  for j in 1:N, i in 1:N
    effratios[i,j] = effs[i]/effs[j]
  end
  bestidx = find((y)->y==maximum(effs),effs)[1]; winner = names[bestidx]
  return Shootout(setups,times,errors,effs,effratios,solutions,names,N,bestidx,winner)
end

function ode_shootoutset(args...;kwargs...)
  warn("ode_shootoutset is deprecated. Use ShootoutSet instead")
  ShootoutSet(args...;kwargs...)
end

function ShootoutSet(probs,setups;probaux=nothing,numruns=20,
                     names=nothing,print_names=false,kwargs...)
  N = length(probs)
  shootouts = Vector{Shootout}(N)
  winners = Vector{String}(N)
  if names == nothing
    names = [string(typeof(setups[i][:alg])) for i=1:length(setups)]
  end
  if probaux == nothing
    probaux = Vector{Dict{Symbol,Any}}(N)
    for i in 1:N
      probaux[i] = Dict{Symbol,Any}()
    end
  end
  for i in eachindex(probs)
    print_names && println(names[i])
    shootouts[i] = Shootout(probs[i],setups;numruns=numruns,names=names,kwargs...,probaux[i]...)
    winners[i] = shootouts[i].winner
  end
  return ShootoutSet(shootouts,probs,probaux,N,winners)
end

Base.length(shoot::Shootout) = shoot.N
Base.size(shoot::Shootout) = length(shoot)
Base.endof(shoot::Shootout) = length(shoot)
Base.getindex(shoot::Shootout,i::Int) = shoot.effs[i]
Base.getindex(shoot::Shootout,::Colon) = shoot.effs

function Base.show(io::IO, shoot::Shootout)
  println(io,"Winner: $(shoot.winner)")
  println(io,"EffRatios: $(shoot.effratios[shoot.bestidx,:])")
end

Base.length(set::ShootoutSet) = set.N
Base.size(set::ShootoutSet) = length(set)
Base.endof(set::ShootoutSet) = length(set)
Base.getindex(set::ShootoutSet,i::Int) = set.shootouts[i]
Base.getindex(set::ShootoutSet,::Colon) = set.shootouts
Base.show(io::IO, set::ShootoutSet) = print(io,"ShootoutSet of $(set.N) shootouts ")


## WorkPrecisions

type WorkPrecision
  prob
  abstols
  reltols
  errors
  times
  name
  N::Int
end

type WorkPrecisionSet
  wps::Vector{WorkPrecision}
  N::Int
  abstols
  reltols
  prob
  setups
  names
end

function WorkPrecision(prob,alg,abstols,reltols,dts=nothing;
                       name=nothing,numruns=20,
                       appxsol=nothing,error_estimate=:final,kwargs...)
  N = length(abstols)
  errors = Vector{Float64}(N)
  times = Vector{Float64}(N)
  if name == nothing
    name = "WP-Alg"
  end
  timeseries_errors = error_estimate ∈ TIMESERIES_ERRORS
  dense_errors = error_estimate ∈ DENSE_ERRORS
  for i in 1:N
    t = 0.0
    if dts == nothing
      sol = solve(prob,alg;kwargs...,abstol=abstols[i],
      reltol=reltols[i],timeseries_errors=timeseries_errors,
      dense_errors = dense_errors) # Compile and get result
      sol = solve(prob,alg,sol.u,sol.t,sol.k;kwargs...,abstol=abstols[i],
      reltol=reltols[i],timeseries_errors=timeseries_errors,
      dense_errors = dense_errors) # Compile and get result
      gc()
    else
      sol = solve(prob,alg;kwargs...,abstol=abstols[i],
      reltol=reltols[i],dt=dts[i],timeseries_errors=timeseries_errors,
      dense_errors = dense_errors) # Compile and get result
      sol = solve(prob,alg,sol.u,sol.t,sol.k;kwargs...,abstol=abstols[i],
      reltol=reltols[i],dt=dts[i],timeseries_errors=timeseries_errors,
      dense_errors = dense_errors) # Compile and get result
      gc()
    end
    t = @elapsed for j in 1:numruns
      if dts == nothing
        solve(prob,alg,sol.u,sol.t,sol.k;kwargs...,
                                  abstol=abstols[i],
                                  reltol=reltols[i],
                                  timeseries_errors=false,
                                  dense_errors = false)
      else
         solve(prob,alg,sol.u,sol.t,sol.k;
                                  kwargs...,abstol=abstols[i],
                                  reltol=reltols[i],dt=dts[i],
                                  timeseries_errors=false,
                                  dense_errors = false)
      end

    end
    t = t/numruns

    if appxsol != nothing
      errsol = calculate_errsol(prob,sol,appxsol)
      errors[i] = mean(errsol.errors[error_estimate])
    else
      errors[i] = mean(sol.errors[error_estimate])
    end
    times[i] = t
  end
  return WorkPrecision(prob,abstols,reltols,errors,times,name,N)
end

# This will only do strong errors
function WorkPrecision(prob::Union{AbstractRODEProblem,AbstractSDEProblem},
                       alg,abstols,reltols,dts=nothing;
                       name=nothing,numruns=20,
                       appxsol=nothing,error_estimate=:final,kwargs...)
  N = length(abstols)
  errors = Vector{Float64}(N)
  times = Vector{Float64}(N)
  local_errors = Vector{Float64}(numruns)
  if name == nothing
    name = "WP-Alg"
  end
  timeseries_errors = error_estimate ∈ TIMESERIES_ERRORS
  dense_errors = error_estimate ∈ DENSE_ERRORS
  for i in 1:N
    t = 0.0
    if dts == nothing
      sol = solve(prob,alg;kwargs...,abstol=abstols[i],
      reltol=reltols[i],timeseries_errors=timeseries_errors,
      dense_errors = dense_errors) # Compile and get result
      sol = solve(prob,alg,sol.u,sol.t;kwargs...,abstol=abstols[i],
      reltol=reltols[i],timeseries_errors=timeseries_errors,
      dense_errors = dense_errors) # Compile and get result
      gc()
    else
      sol = solve(prob,alg;kwargs...,abstol=abstols[i],
      reltol=reltols[i],dt=dts[i],timeseries_errors=timeseries_errors,
      dense_errors = dense_errors) # Compile and get result
      sol = solve(prob,alg,sol.u,sol.t;kwargs...,abstol=abstols[i],
      reltol=reltols[i],dt=dts[i],timeseries_errors=timeseries_errors,
      dense_errors = dense_errors) # Compile and get result
      gc()
    end
    t = @elapsed for j in 1:numruns
      if dts == nothing
        solve(prob,alg,sol.u,sol.t;kwargs...,
                                  abstol=abstols[i],
                                  reltol=reltols[i],
                                  timeseries_errors=false,
                                  dense_errors = false)
      else
        solve(prob,alg,sol.u,sol.t;
                                  kwargs...,abstol=abstols[i],
                                  reltol=reltols[i],dt=dts[i],
                                  timeseries_errors=false,
                                  dense_errors = false)
      end
      if appxsol != nothing
        errsol = calculate_errsol(prob,sol,appxsol)
        local_errors[j] = errsol.errors[error_estimate]
      else
        local_errors[j] = sol.errors[error_estimate]
      end
    end
    t = t/numruns

    errors[i] = mean(local_errors)
    times[i] = t
  end
  return WorkPrecision(prob,abstols,reltols,errors,times,name,N)
end

function calculate_errsol(prob,sol::AbstractODESolution,appxsol_setup::Dict)
  true_sol = solve(prob,appxsol_setup[i][:alg];appxsol_setup[i]...)
  appxtrue(sol,true_sol)
end

function calculate_errsol(prob::AbstractSDEProblem,sol::AbstractRODESolution,appxsol_setup::Dict)
  prob2 = SDEProblem(prob.f,prob.g,prob.u0,prob.tspan,noise=NoiseWrapper(sol.W))
  true_sol = solve(prob2,appxsol_setup[i][:alg];appxsol_setup[i]...)
  appxtrue(sol,true_sol)
end

function calculate_errsol(prob,sol::AbstractODESolution,true_sol::AbstractTimeseriesSolution)
  appxtrue(sol,true_sol)
end

function calculate_errsol(prob::MonteCarloProblem,sol,true_sol)
  appxtrue(sol,true_sol)
end

function WorkPrecisionSet(prob,abstols,reltols,setups;numruns=20,
                          print_names=false,names=nothing,appxsol=nothing,kwargs...)
  N = length(setups)
  wps = Vector{WorkPrecision}(N)
  if names == nothing
    names = [string(typeof(setups[i][:alg])) for i=1:length(setups)]
  end
  for i in 1:N
    print_names && println(names[i])
    if haskey(setups[i],:dts)
      wps[i] = WorkPrecision(prob,setups[i][:alg],abstols,reltols,setups[i][:dts];
                                 numruns=numruns,appxsol=appxsol,
                                 name=names[i],kwargs...,setups[i]...)
    else
      wps[i] = WorkPrecision(prob,setups[i][:alg],abstols,reltols;
                                 numruns=numruns,appxsol=appxsol,
                                 name=names[i],kwargs...,setups[i]...)
    end
  end
  return WorkPrecisionSet(wps,N,abstols,reltols,prob,setups,names)
end

Base.length(wp::WorkPrecision) = wp.N
Base.size(wp::WorkPrecision) = length(wp)
Base.endof(wp::WorkPrecision) = length(wp)
Base.getindex(wp::WorkPrecision,i::Int) = wp.times[i]
Base.getindex(wp::WorkPrecision,::Colon) = wp.times

function Base.show(io::IO, wp::WorkPrecision)
  println(io,"Name: $(wp.name)")
  println(io,"Times: $(wp.times)")
  println(io,"Errors: $(wp.errors)")
end

Base.length(wp_set::WorkPrecisionSet) = wp_set.N
Base.size(wp_set::WorkPrecisionSet) = length(wp_set)
Base.endof(wp_set::WorkPrecisionSet) = length(wp_set)
Base.getindex(wp_set::WorkPrecisionSet,i::Int) = wp_set.wps[i]
Base.getindex(wp_set::WorkPrecisionSet,::Colon) = wp_set.wps

Base.show(io::IO, wp_set::WorkPrecisionSet) = println(io,"WorkPrecisionSet of $(wp_set.N) wps")
