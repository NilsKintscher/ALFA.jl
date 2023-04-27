
@recipe function f(
    L::Lattice{1,T};
    xmin = -3,
    xmax = 3,
    draw_basis = true,
) where {T}
    seriescolor --> :gray
    linewidth = get(plotattributes, :linewidth, 0.5)

    ylims --> (-.5, 0.5)

    @series begin
        x = [L.A * xmin; L.A * xmax]
        y = zeros(size(x))
        linewidth := linewidth
        label := ""
        x, y
    end
    @series begin
        #xy = hcat([(L.A * [i])' for i = xmin:xmax]...)
        seriestype --> :scatter

        markershape --> :vline
        markeralpha := 1
        markerstrokealpha := 1
        markersize --> 10
        markerstrokewidth --> 5
        markerstrokecolor --> :black

        x = hcat([L.A * i for i = xmin:xmax]...)
        y = zeros(size(x))
        label := ""
        primary := false
        x, y
    end
    if draw_basis
        seriescolor := :black
        arrow := true
        label --> ""
        linewidth := 2 * linewidth
        @series begin
            x = [0; L.A[1, 1]]
            y = [0; 0]
            x, y
        end
    end
end

@recipe function f(
    L::Lattice{2,T};
    xmin = -3,
    xmax = 3,
    draw_basis = true,
) where {T}
    seriescolor --> :gray
    linewidth = get(plotattributes, :linewidth, 0.5)
    @series begin
        xy = hcat([(L.A * [i i; xmin * 1.1 xmax * 1.1])' for i = xmin:xmax]...)
        x = xy[:, 1:2:end]
        y = xy[:, 2:2:end]
        label := ""
        primary := false
        linewidth := linewidth
        x, y
    end
    @series begin
        xy = hcat([(L.A * [xmin * 1.1 xmax * 1.1; i i])' for i = xmin:xmax]...)
        x = xy[:, 1:2:end]
        y = xy[:, 2:2:end]
        linewidth := linewidth
        label := ""
        primary := false #
        x, y
    end
    if draw_basis
        seriescolor := :black
        arrow := true
        label --> ""
        linewidth := 2 * linewidth
        @series begin
            x = [0 0; L.A[1, :]']
            y = [0 0; L.A[2, :]']
            x, y
        end
    end
end

@recipe function f(
    C::Crystal;
    xmin = -3,
    xmax = 3,
    draw_basis = true,
    plot_domain = true,
    plot_codomain = true,
    plot_lattice = true,
)
    if C.n == 1
        ylims := (-.5, 0.5)
    end
    if plot_lattice
        @series begin
            C.L
        end
    end
    pos_fractional = Iterators.product(Iterators.repeated(xmin:xmax, C.dim)...) # fractional positions

    if plot_codomain
        @series begin  #codomain
            markercolor --> :white
            label := ""
            markershape --> :pentagon
            markeralpha := 1
            markerstrokealpha := 1
            markersize --> 10
            markerstrokewidth --> 1
            markerstrokecolor --> :black
            plot_domain := false
            C, pos_fractional#, domain=false
        end
    end
    if plot_domain
        @series begin #domain
            markercolor --> :salmon
            markerstrokealpha := 0
            label := ""
            markershape --> :diamond
            plot_domain := true
            C, pos_fractional#, domain=true
        end
    end

end

@recipe function f(
    C::Crystal{1,T},
    pos_fractional;
    plot_domain = true,
) where {T}

    seriestype --> :scatter
    markercolor --> :orange
    label := ""
    markershape --> :dtriangle
    if plot_domain
        s = C.Domain
    else
        s = C.Codomain
    end

    xy =  C.L.A * collect(Iterators.flatten(pos_fractional))'
    for pos in s # eachslice(s, dims = 1)
        @series begin
            xy[1, :] .+ pos[1], zeros(size(xy[1, :]))
        end
    end
end

@recipe function f(
    C::Crystal{2,T},
    pos_fractional;
    plot_domain = true,
) where {T}
    seriestype --> :scatter
    markercolor --> :orange
    label := ""
    markershape --> :dtriangle
    if plot_domain
        s = C.Domain
    else
        s = C.Codomain
    end

    xy = C.L.A * reshape(collect(Iterators.flatten(pos_fractional)), 2, :)
    for pos in s # eachslice(s, dims = 1)
        @series begin
            xy[1, :] .+ pos[1], xy[2, :] .+ pos[2]
        end
    end
end





@recipe function f(S::CrystalOperator{1,T}; threshhold = 1e-12) where {T}

    ylims := (-.5, 0.5)
    m = collect(S.M)
    x = vcat(transpose([x.pos for x in S.M])...)
    # extrema of positions
    (xmin, xmax) = extrema(x)
    # extrema of matrix entries
    mmin = min([min(real(mu.mat)...) for mu in m]...)
    mmax = max([max(real(mu.mat)...) for mu in m]...)

    if mmin ≈ mmax

        mmax = mmax + abs(0.1 * mmax)
        mmin = mmin - abs(0.1 * mmin)
    end

    @series begin # plot lattice
        xmin --> xmin
        xmax --> xmax + 1
        S.C.L
    end


    @series begin #
        plot_lattice := false
        plot_domain := false
        plot_codomain := true

        xmin --> 0
        xmax --> 0
        S.C # codomain
    end
    @series begin
        plot_lattice := false
        plot_domain := true
        plot_codomain := false
        xmin --> xmin
        xmax --> xmax
        S.C
    end

    ### colorbar
    @series begin
        label --> ""
        seriescolor := :viridis
        line_z := range(mmin, stop = mmax, length = 10)
        [NaN], [NaN]
    end

    for m in S.M #  pos != 0.
        coord_cartesian = S.C.L.A * m.pos
        for (it_d, sd) in enumerate(S.C.Domain) #eachslice(S.C.Domain, dims = 1))
            for (it_c, sc) in enumerate(S.C.Codomain) #eachslice(S.C.Codomain, dims = 1))
                val = real(m.mat[it_c, it_d])
                if abs(val) > threshhold
                    p0 = [coord_cartesian + sd..., 0]
                    p2 = [sc..., 0]
                    if norm(p0 - p2) ≈ 0
                    else
                        linecolor := get(
                            ColorSchemes.viridis,
                            Float64(-(mmin - val) / (mmax - mmin)),
                        ) # linear interpolation
                        label := ""
                        p2 = [sc..., 0]
                        d = -p0 + p2
                        mid = p0 .+ d ./ 2
                        nd = [-d[2], d[1]]

                        p1 = mid + 0.3 * [-d[2], d[1]]#/norm(d)*r

                        B(t) = (1 - t)^2 * p0 + 2 * (1 - t) * t * p1 + t^2 * p2
                        vv = vcat([
                            B(t) for t in range(0, stop = 0.95, step = 0.01)
                        ]'...)

                        vv = Array(vv)
                        @series begin
                            arrow := true
                            linewidth := 1.75
                            vv[:, 1], vv[:, 2]
                        end
                    end
                end
            end
        end
    end

    for m in S.M #  pos != 0.
        coord_cartesian = S.C.L.A * m.pos
        for (it_d, sd) in enumerate(S.C.Domain) #eachslice(S.C.Domain, dims = 1))
            for (it_c, sc) in enumerate(S.C.Codomain) #eachslice(S.C.Codomain, dims = 1))
                val = real(m.mat[it_c, it_d])
                #if abs(val) > threshhold
                p0 = [coord_cartesian + sd, 0]
                p2 = [sc, 0]
                if norm(p0 - p2) ≈ 0
                    @series begin
                        seriestype := :scatter
                        markercolor := get(
                            ColorSchemes.viridis,
                            Float64(-(mmin - val) / (mmax - mmin)),
                        ) # linear interpolation
                        label := ""
                        markershape --> :pentagon
                        markeralpha := 1
                        markersize --> 5 # 10
                        markerstrokewidth --> 0
                        #marker_z = -4:-4
                        pos = coord_cartesian + sd
                        [pos[1]], [0]
                    end
                end
                #end
            end
        end
    end
end


@recipe function f(S::CrystalOperator{2,T}; threshhold = 1e-12) where {T}
    m = collect(S.M)
    x = vcat(transpose([x.pos for x in S.M])...)
    # extrema of positions
    (xmin, xmax) = extrema(x)
    # extrema of matrix entries
    mmin = min([min(real(mu.mat)...) for mu in m]...)
    mmax = max([max(real(mu.mat)...) for mu in m]...)

    if mmin ≈ mmax

        mmax = mmax + abs(0.1 * mmax)
        mmin = mmin - abs(0.1 * mmin)
    end

    @series begin # plot lattice
        xmin --> xmin
        xmax --> xmax + 1
        S.C.L
    end


    @series begin #
        plot_lattice := false
        plot_domain := false
        plot_codomain := true

        xmin --> 0
        xmax --> 0
        S.C # codomain
    end
    @series begin
        plot_lattice := false
        plot_domain := true
        plot_codomain := false
        xmin --> xmin
        xmax --> xmax
        S.C
    end

    ### colorbar
    @series begin
        label --> ""
        seriescolor := :viridis
        line_z := range(mmin, stop = mmax, length = 10)
        [NaN], [NaN]
    end


    for m in S.M #  pos != 0.
        coord_cartesian = S.C.L.A * m.pos
        for (it_d, sd) in enumerate(S.C.Domain) #eachslice(S.C.Domain, dims = 1))
            for (it_c, sc) in enumerate(S.C.Codomain) #eachslice(S.C.Codomain, dims = 1))
                val = real(m.mat[it_c, it_d])
                if abs(val) > threshhold
                    p0 = coord_cartesian + sd
                    p2 = sc
                    if norm(p0 - p2) ≈ 0
                    else
                        linecolor := get(
                            ColorSchemes.viridis,
                            Float64(-(mmin - val) / (mmax - mmin)),
                        ) # linear interpolation
                        label := ""
                        p2 = sc
                        d = -p0 + p2
                        mid = p0 .+ d ./ 2
                        nd = [-d[2], d[1]]
                        p1 = mid + 0.3 * [-d[2], d[1]]#/norm(d)*r
                        B(t) = (1 - t)^2 * p0 + 2 * (1 - t) * t * p1 + t^2 * p2
                        vv = vcat([
                            B(t) for t in range(0, stop = 0.9, step = 0.01)
                        ]'...)

                        vv = Array(vv)
                        @series begin
                            arrow := true
                            linewidth := 1.75
                            vv[:, 1], vv[:, 2]
                        end
                    end
                end
            end
        end
    end

    for m in S.M #  pos != 0.
        coord_cartesian = S.C.L.A * m.pos
        for (it_d, sd) in enumerate(S.C.Domain) #eachslice(S.C.Domain, dims = 1))
            for (it_c, sc) in enumerate(S.C.Codomain) #eachslice(S.C.Codomain, dims = 1))
                val = real(m.mat[it_c, it_d])
                #if abs(val) > threshhold
                p0 = coord_cartesian + sd
                p2 = sc
                if norm(p0 - p2) ≈ 0
                    @series begin
                        seriestype := :scatter
                        markercolor := get(
                            ColorSchemes.viridis,
                            Float64(-(mmin - val) / (mmax - mmin)),
                        ) # linear interpolation
                        label := ""
                        markershape --> :pentagon
                        markeralpha := 1
                        markersize --> 5 # 10
                        markerstrokewidth --> 0
                        #marker_z = -4:-4
                        pos = coord_cartesian + sd
                        [pos[1]], [pos[2]]
                    end
                end
                #end
            end
        end
    end
end


@userplot plotSpectrum

@recipe function f(h::plotSpectrum; N = 20, zfilter = nothing, modifier = abs, surfaceIndexFirst = lastindex, surfaceIndexSecond = lastindex )

    @assert h.args[1] isa CrystalOperator || h.args[1] isa OperatorComposition "input must be a CrystalOperator or OperatorComposition"
    S = h.args[1]

    if S.C.n == 1
        x = range(0, stop = 1, length = N + 1)[1:end-1]
        maxval = -Inf
        maxx = NaN
        
        function f1d(x, surfaceIdx = surfaceIdx)
            dAxy = S.C.L.dA * [x]
            try
                evals = sort(modifier.(ALFA.eigvals(S, dAxy)))

                z =  evals[surfaceIdx]
            catch
                z = NaN
            end
            if zfilter !== nothing
                if z < zfilter[1] || z > zfilter[2]
                    z = NaN
                end
            end
            if z > maxval
                maxval = z
                maxx = x
            end
            return z
        end

          
        num_evals = 1:S.C.size_codomain

        if typeof(surfaceIndexFirst) <: Int 
            idxRangeStart = surfaceIndexFirst
        else
            idxRangeStart = surfaceIndexFirst(num_evals)
        end
        if typeof(surfaceIndexSecond) <: Int 
            idxRangeEnd = surfaceIndexSecond
        else
            idxRangeEnd = surfaceIndexSecond(num_evals)
        end

        for idx in idxRangeStart:idxRangeEnd


        zv = [f1d(xx, idx) for xx in x]

        seriescolor --> :viridis
        @series begin
            #seriestype := :surface
            label := ""
            title := ""
            x, zv

        end
        @series begin
            seriestype := :scatter
            markercolor := :red
            label -> "max(z) = " * fmt(FormatSpec(".3e"), maxval)
            [maxx], [maxval]
        end
    end 
    elseif S.C.n == 2
        x = y = range(0, stop = 1, length = N + 1)[1:end-1]
        maxval = -Inf
        maxx = NaN
        maxy = NaN
        function f2d(x, y, surfaceIdx = surfaceIdx)
            dAxy = S.C.L.dA * [x, y]
            try
                evals = sort(modifier.(ALFA.eigvals(S, dAxy)))

                z =  evals[surfaceIdx]
            catch
                z = NaN
            end 
            if zfilter !== nothing
                if z < zfilter[1] || z > zfilter[2]
                    z = NaN
                end
            end
            if z > maxval
                maxval = z
                maxx = x
                maxy = y
            end
            return z
        end
        
        num_evals = 1:S.C.size_codomain

        if typeof(surfaceIndexFirst) <: Int 
            idxRangeStart = surfaceIndexFirst
        else
            idxRangeStart = surfaceIndexFirst(num_evals)
        end
        if typeof(surfaceIndexSecond) <: Int 
            idxRangeEnd = surfaceIndexSecond
        else
            idxRangeEnd = surfaceIndexSecond(num_evals)
        end

        for idx in idxRangeStart:idxRangeEnd

            zv = [f2d(xx, yy, idx) for xx in x for yy in y]

            layout := (1, 2)
            seriescolor --> :viridis
            @series begin
                subplot := 1
                seriestype := :surface
                title := ""
                x, y, zv

            end

            @series begin
                subplot := 2
                seriestype := :contourf
                x, y, zv
            end

            @series begin
                subplot := 2
                seriestype := :scatter
                markercolor := :red
                markersize := 4
                label := "max(z) = " * fmt(FormatSpec(".3e"), maxval)
                [maxx], [maxy]
            end
     end

    end
end


"""
plotspectrum(L, N = 20, zfilter = nothing, modifier = abs, surfaceIndexFirst = lastindex, surfaceIndexSecond = lastindex )

Plots the spectrum of L along its dual lattice with N points in each direction.
- Implemented only for 1d and 2d.

For each dual lattice point you may have multiple complex eigenvalues. In order to plot them, you need to modify them by taking for example the absolute part of the eigenvalues.
Furthermore, if you have multiple eigenvalues per dual lattice point, you are able to control what surfaces you want to plot via the arguuments surfaceIndexFirst and surfaceIndexSecond.

You can use zfilter in order to show only the part of the plot between zfilter[1] and zfilter[2]
# Example
```jldoctest
julia> using ALFA
julia> using Plots

julia> L = ALFA.gallery.graphene_tight_binding()
Lattice Basis: ALFA.Lattice{2, Float64}([1.5 1.5; 0.8660254037844386 -0.8660254037844386])
Domain: 2-element Vector{StaticArraysCore.SVector{2, Float64}}:
 [1.0, 0.0]
 [2.0, 0.0]
Codomain: 2-element Vector{StaticArraysCore.SVector{2, Float64}}:
 [1.0, 0.0]
 [2.0, 0.0]
Multiplier: 9-element Vector{ALFA.Multiplier}:
 ALFA.Multiplier{2}([-1, -1], [0 0; 0 0])
 ALFA.Multiplier{2}([-1, 0], [0 -1; 0 0])
 ALFA.Multiplier{2}([-1, 1], [0 0; 0 0])
 ALFA.Multiplier{2}([0, -1], [0 -1; 0 0])
 ALFA.Multiplier{2}([0, 0], [0 -1; -1 0])
 ALFA.Multiplier{2}([0, 1], [0 0; -1 0])
 ALFA.Multiplier{2}([1, -1], [0 0; 0 0])
 ALFA.Multiplier{2}([1, 0], [0 0; -1 0])
 ALFA.Multiplier{2}([1, 1], [0 0; 0 0])

# we have two eigenvalues per dual lattice point.
julia> plotspectrum(L) # plots the absolute part of the largest absolute eigenvalue

julia> plotspectrum(L, modifier=real, surfaceIndexFirst=firstindex, surfaceIndexSecond=lastindex) # plots the real part of all eigenvalues
julia> plotspectrum(L, modifier=real, surfaceIndexFirst=1, surfaceIndexSecond=2) # equivalent to the functioncall above

julia> plotspectrum(L, modifier=real, surfaceIndexFirst=firstindex, surfaceIndexSecond=firstindex) # plots the real part of the smallest eigenvalues

julia> plotspectrum(L, zfilter=[2, 3], modifier=abs, surfaceIndexFirst=lastindex, surfaceIndexSecond=lastindex) # plots the abs part of largest eigenvalues that lie between 2 and 3.

julia> plotspectrum(L, zfilter=[-1,1], modifier=imag, surfaceIndexFirst=1, surfaceIndexSecond=2, zlim=[-1,1]) # modifying z axis of the plot is also possible.

"""
plotspectrum