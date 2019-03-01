@recipe function f(sim::ConvergenceSimulation)
  if any((x)->ndims(x)>1,values(sim.errors)) #Monte Carlo
    vals = [typeof(x)<:AbstractVector ? x : vec(mean(x,1)) for x in values(sim.errors)]
  else #Deterministic
    vals = [x for x in values(sim.errors)]
  end
  seriestype --> :path
  label  --> reshape([string(key) for key in keys(sim.errors)],1,length(keys(sim.errors)))
  xguide  --> "Convergence Axis"
  yguide  --> "Error"
  linewidth --> 3
  xscale --> :log10
  yscale --> :log10
  marker --> :auto
  sim.convergence_axis, vals
end

@recipe function f(shoot::Shootout)
  seriestype --> :bar
  legend := false
  xguide --> "Algorithms"
  yguide --> "Efficiency"
  shoot.names,shoot.effs
end

@recipe function f(wp::WorkPrecision)
  seriestype --> :path
  label -->  wp.name
  linewidth --> 3
  yguide --> "Time (s)"
  xguide --> "Error"
  xscale --> :log10
  yscale --> :log10
  marker --> :auto
  wp.errors,wp.times
end

@recipe function f(wp_set::WorkPrecisionSet;
                   plot_sample_error = wp_set.error_estimate ∈ WEAK_ERRORS)
  seriestype --> :path
  linewidth --> 3
  yguide --> "Time (s)"
  xguide --> "Error"
  xscale --> :log10
  yscale --> :log10
  marker --> :auto
  errors = Vector{Any}(undef,0)
  times = Vector{Any}(undef,0)
  for i in 1:length(wp_set)
    push!(errors,wp_set[i].errors)
    push!(times,wp_set[i].times)
  end
  if plot_sample_error
    @series begin
      linestyle := :dash
      label := "Sample Error: $(wp_set.numruns)"
      color := :red
      xs = [wp_set.sample_error,wp_set.sample_error]
      ys = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
      xs,ys
    end
    label -->  reshape(["Sample Error: $(wp_set.numruns)";wp_set.names],1,1+length(wp_set))
  else
    label -->  reshape(wp_set.names,1,length(wp_set))
  end

  errors,times
end

@recipe function f(tab::ODERKTableau;dx=1/100,dy=1/100,xlim=[-6,1],ylim=[-5,5])
  x = xlim[1]:dx:xlim[2]
  y = ylim[1]:dy:ylim[2]
  f = (u,v)-> abs(stability_region(u+v*im,tab))<1
  seriestype --> :contour
  fill --> true
  cbar --> false
  x,y,f
end
