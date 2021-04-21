@recipe function f(
    CV::CrystalVector{N,T,outerdim, innerdim}
    ) where {N,T,outerdim, innerdim}
    xy = hcat(vcat([[y+CV.CT.C.L.A*x for y in CV.CT.C.Domain] for x in  CV.CT.coords]...)...)
    x = xy[1,:]
    z = vcat(CV.v...)
    if N==2
    y = xy[2,:]
    @series begin
        seriescolor --> :viridis
        seriestype := :scatter #:surface
        linewidth --> 1.0
        #opacity --> 0.5#
        marker_z := z
        x, y#, z
        # x, y, z  # for surface
    end
elseif N==1
    @series begin
        x, z
    end
end
end

@recipe function f(
    CT::CrystalTorus
    )
    @series begin
        xmin := 0
        xmax := 1
        CT.Z
    end
    @series begin
        CT.C, CT.coords
    end
end


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

@recipe function f(h::plotSpectrum; N = 20, zfilter = nothing)

    @assert h.args[1] isa CrystalOperator || h.args[1] isa OperatorComposition "input must be a CrystalOperator or OperatorComposition"
    S = h.args[1]

    if S.C.n == 1
        x = range(0, stop = 1, length = N + 1)[1:end-1]
        maxval = -Inf
        maxx = NaN

        function f1d(x)
            dAxy = S.C.L.dA * [x]
            try
                z = abs(ALFA.eigvals(S, dAxy)[end])
            catch
                z = NaN
            end
            if zfilter != nothing
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

        zv = [f1d(xx) for xx in x]

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
    elseif S.C.n == 2
        x = y = range(0, stop = 1, length = N + 1)[1:end-1]
        maxval = -Inf
        maxx = NaN
        maxy = NaN

        function f2d(x, y)
            dAxy = S.C.L.dA * [x, y]
            try
                z = abs(ALFA.eigvals(S, dAxy)[end])
            catch
                z = NaN
            end
            if zfilter != nothing
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

        zv = [f2d(xx, yy) for xx in x for yy in y]

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
