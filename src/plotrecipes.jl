
@recipe function f(L::Lattice; xmin = -3, xmax = 3, draw_basis = true)
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
        primary := false #
        linewidth := linewidth
        label := ""
        x, y
    end
    if draw_basis
        color := :black
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
            domain := false
            C, pos_fractional#, domain=false
        end
    end
    if plot_domain
        @series begin #domain
            markercolor --> :salmon
            markerstrokealpha := 0
            label := ""
            markershape --> :diamond
            domain := true
            C, pos_fractional#, domain=true
        end
    end

end

@recipe function f(C::Crystal, pos_fractional; domain = true)
    # #TODO:  check this:
    # domain --> true
    # codomain --> !domain
    # #

    seriestype --> :scatter
    markercolor --> :orange
    label := ""
    markershape --> :dtriangle
    if domain
        s = C.Domain
    else
        s = C.Codomain
    end

    xy = C.L.A * reshape(collect(Iterators.flatten(pos_fractional)), 2, :)
    for pos in eachslice(s, dims = 1)
        @series begin
            xy[1, :] .+ pos[1], xy[2, :] .+ pos[2]
        end
    end
end







@recipe function f(S::CrystalOperator; threshhold = 1e-15)
    m = collect(S.M)
    x = vcat(transpose([x.pos for x in S.M])...)
    # extrema of positions
    (xmin, xmax) = extrema(x)
    # extrema of matrix entries
    mmin = min([min(real(mu.mat)...) for mu in m]...)
    mmax = max([max(real(mu.mat)...) for mu in m]...)
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
        c := :viridis
        line_z := mmin:mmax
        [], []
    end

    for m in S.M #  pos != 0.
        coord_cartesian = S.C.L.A * m.pos
        #println("mat: $(m.mat)")
        for (it_d, sd) in enumerate(eachslice(S.C.Domain, dims = 1))
            for (it_c, sc) in enumerate(eachslice(S.C.Codomain, dims = 1))
                val = real(m.mat[it_c, it_d])
                if abs(val) > threshhold
                    p0 = coord_cartesian + sd
                    p2 = sc
                    if norm(p0 - p2) ≈ 0

                    else
                        linecolor := get(
                            ColorSchemes.viridis,
                            -(mmin - val) / (mmax - mmin),
                        ) # linear interpolation
                        label := ""

                        # # draw arrow.
                        #println("coord_cartesian: $coord_cartesian , sd: $sd")
                        #println("p0")
                        #println(p0)
                        p2 = sc
                        #println("p2")
                        #println(p2)
                        d = -p0 + p2
                        #println("d")
                        #println(d)

                        mid = p0 .+ d ./ 2
                        #println("mid")
                        #println(mid)
                        nd = [-d[2], d[1]]
                        p1 = mid + 0.3 * [-d[2], d[1]]#/norm(d)*r
                        #println("p1")
                        #println(p1)
                        B(t) = (1 - t)^2 * p0 + 2 * (1 - t) * t * p1 + t^2 * p2
                        ## shorten by 10pt

                        #
                        vv = vcat([
                            B(t) for t in range(0, stop = 1, step = 0.01)
                        ]'...)
                        @series begin
                            #arrow := true
                            linewidth := 2
                            vv[:, 1], vv[:, 2]
                        end
                    end
                end
            end
        end

    end


    for m in S.M #  pos != 0.
        coord_cartesian = S.C.L.A * m.pos
        for (it_d, sd) in enumerate(eachslice(S.C.Domain, dims = 1))
            for (it_c, sc) in enumerate(eachslice(S.C.Codomain, dims = 1))
                val = real(m.mat[it_d, it_c])
                #if abs(val) > threshhold
                p0 = coord_cartesian + sd
                p2 = sc
                if norm(p0 - p2) ≈ 0
                    println(val)
                    println(mmin)
                    println(mmax)
                    print(-(mmin - val) / (mmax - mmin))
                    @series begin
                        seriestype := :scatter
                        markercolor := get(
                            ColorSchemes.viridis,
                            -(mmin - val) / (mmax - mmin),
                        ) # linear interpolation
                        #markercolor := get(ColorSchemes.viridis,0.5) # linear interpolation
                        #clim := -5:-4
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
