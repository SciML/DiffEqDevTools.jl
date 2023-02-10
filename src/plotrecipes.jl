@recipe function f(sim::ConvergenceSimulation)
    if any((x) -> ndims(x) > 1, values(sim.errors)) #Monte Carlo
        vals = [typeof(x) <: AbstractVector ? x : vec(mean(x, 1))
                for x in values(sim.errors)]
    else #Deterministic
        vals = [x for x in values(sim.errors)]
    end
    seriestype --> :path
    label -->
    reshape([string(key) for key in keys(sim.errors)], 1, length(keys(sim.errors)))
    xguide --> "Convergence Axis"
    yguide --> "Error"
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
    shoot.names, shoot.effs
end

@recipe function f(wp::WorkPrecision)
    seriestype --> :path
    label --> wp.name
    linewidth --> 3
    yguide --> "Time (s)"
    xguide --> "Error"
    xscale --> :log10
    yscale --> :log10
    markershape --> :auto
    wp.errors, wp.times
end

@recipe function f(wp_set::WorkPrecisionSet; view = :benchmark, color = nothing)
    if view == :benchmark
        seriestype --> :path
        linewidth --> 3
        yguide --> "Time (s)"
        xguide --> "Error"
        xscale --> :log10
        yscale --> :log10
        markershape --> :auto
        errors = Vector{Any}(undef, 0)
        times = Vector{Any}(undef, 0)
        for i in 1:length(wp_set)
            push!(errors, wp_set[i].errors)
            push!(times, wp_set[i].times)
        end
        label --> reshape(wp_set.names, 1, length(wp_set))
        return errors, times
    elseif view == :dt_convergence
        idts = filter(i -> haskey(wp_set.setups[i], :dts), 1:length(wp_set))
        length(idts) > 0 ||
            throw(ArgumentError("Convergence with respect to Δt requires runs with fixed time steps"))
        dts = Vector{Any}(undef, 0)
        errors = Vector{Any}(undef, 0)
        ps = Vector{Any}(undef, 0)
        convs = Vector{Any}(undef, 0)
        for i in idts
            push!(dts, wp_set.setups[i][:dts])
            push!(errors, wp_set[i].errors)
            lc, p = [one.(dts[end]) log.(dts[end])] \ log.(errors[end])
            push!(ps, p)
            push!(convs, exp(lc) * dts[end] .^ p)
        end
        names = wp_set.names[idts] .*
                map(p -> " (Δtᵖ order p=$(round(p, sigdigits=2)))", ps)
        if color === nothing
            color = reshape(1:length(idts), 1, :)
        end
        @series begin
            seriestype --> :path
            linestyle --> :dash
            linewidth --> 1
            color --> color
            xscale --> :log10
            yscale --> :log10
            markersize --> 0
            label --> nothing
            return dts, convs
        end
        @series begin
            seriestype --> :path
            linewidth --> 3
            color --> color
            yguide --> "Error"
            xguide --> "Δt"
            xscale --> :log10
            yscale --> :log10
            markershape --> :auto
            label --> reshape(names, 1, length(idts))
            return dts, errors
        end
    else
        throw(ArgumentError("view argument `$view` not implemented"))
    end
end

@recipe function f(tab::ODERKTableau; dx = 1 / 100, dy = 1 / 100, order_star = false,
                   embedded = false)
    xlims = get(plotattributes, :xlims, (-6, 1))
    ylims = get(plotattributes, :ylims, (-5, 5))
    x = xlims[1]:dx:xlims[2]
    y = ylims[1]:dy:ylims[2]

    if order_star
        f = (u, v) -> abs(stability_region(u + v * im, tab; embedded = embedded) /
                          exp(u + v * im)) < 1
    else
        f = (u, v) -> abs(stability_region(u + v * im, tab; embedded = embedded)) < 1
    end
    seriestype --> :contour
    fill --> true
    colorbar --> false
    seriescolor --> :grays
    aspect_ratio --> 1
    x, y, f
end
