struct CrystalOperator
    C::Crystal ## dimension L.dim dictates size of Multiplier.pos and dim of structure elements.
    M::SortedSet{Multiplier}# Array{Multiplier,1}
    _CompatibilityCheckOnly::Bool
    function CrystalOperator(
        C::Crystal,
        M::SortedSet{Multiplier},
        _CompatibilityCheckOnly::Bool,
    )
        _sanitycheck(C, M)
        new(C, M, _CompatibilityCheckOnly)
    end
end

function CrystalOperator(
    C = nothing,
    M = nothing,
    _CompatibilityCheckOnly = false,
)
    if C == nothing
        C = Crystal()
    end
    if M == nothing
        M = SortedSet{Multiplier}()
    end
    return CrystalOperator(C, M, _CompatibilityCheckOnly)
end

function _sanitycheck(C::Crystal, m::Multiplier) # check dimensionality and size:
    #dim / pos:
    @assert m.dim == C.dim "Multiplier dimension m.dim=length(m.pos)=$(m.dim) must be equal to Lattice dimension S.C.dim=size(S.C.L,1)=$(C.dim)"
    #size domain, codomain
    @assert m.size_domain == C.size_domain "Multiplier domain size m.size_domain=size(m.mat,2)=$(m.size_domain) must be equal to size of the (structure element of the) Domain C.size_domain=size(C.Domain,1)=$(C.size_domain)"
    @assert m.size_codomain == C.size_codomain "Multiplier codomain size m.size_codomain=size(m.mat,1)=$(m.size_codomain) must be equal to size of the (structure element of the) Domain C.size_codomain=size(C.Codomain,1)=$(C.size_codomain)"
    return true
end

function _sanitycheck(C::Crystal, M::SortedSet{Multiplier})
    for m in M
        _sanitycheck(C, m)
    end
end

function Base.push!(S::CrystalOperator, m::Multiplier)
    _sanitycheck(S.C, m)
    push!(S.M, m)
    return S
end


function Base.show(io::IO, mime::MIME"text/plain", s::SortedSet{Multiplier})
    show(io, mime, collect(s))
end

function Base.show(io::IO, mime::MIME"text/plain", o::CrystalOperator)
    show(io, mime, o.C)
    print(io, "\nMultiplier: ")
    show(io, mime, o.M)
end

function normalize(S::CrystalOperator)
    (dn, ds, dp) = ShiftIntoUnitCell(S.C.Domain, S.C.L)
    (cn, cs, cp) = ShiftIntoUnitCell(S.C.Codomain, S.C.L)

    # ...

end

function wrtLattice(S::CrystalOperator, A::Matrix)
    #### TODO: Test this function, when Plot function exists. Need to construct test cases with alfa.py.
    t = ElementsInQuotientSpace(S.C.A, A, fractional = true)
    newDomain = vcat(transpose([x + y for x in eachslice(t, dims = 1) for y in eachslice(
        S.C.Domain,
        dims = 1,
    )])...)
    newCodomain = vcat(transpose([x + y for x in eachslice(t, dims = 1) for y in eachslice(
        S.C.Codomain,
        dims = 1,
    )])...)

    m_old = collect(S.M)
    y_old = vcat(transpose([x.pos for x in S.M])...)
    Ay_old = transpose(S.C.L.A * transpose(y_old))
    # find all unique combinations of floor(A\S.C.L.A*(y_old[i] + t[j] - t[k]))
    SS = SortedSet{Array{Int,1}}()
    for j in Iterators.product(
        eachslice(y_old, dims = 1),
        eachslice(t, dims = 1),
        eachslice(t, dims = 1),
    )
        y = A \ (S.C.L.A * (j[1] + j[2] - j[3]))
        map!(
            x -> isapprox(x, round(x), rtol = alfa_rtol, atol = alfa_atol) ?
                 round(x) : floor(x),
            y,
            y,
        )
        push!(SS, y)
    end
    y_new = vcat(transpose(collect(SS))...)

    Ay_new = transpose(A * transpose(y_new)) # convert to cartesian coordinate of original lattice.

    # get coordinate of new lattice
    print(Ay_new)
    # assign multipliers.

    brs = S.C.size_codomain #block row size
    bcs = S.C.size_domain #block column size

    Cnew = Crystal(A, newDomain, newCodomain)
    op = CrystalOperator(Cnew)
    for (it_y, y) in enumerate(eachslice(Ay_new, dims = 1))
        mm = zeros(
            typeof(first(S.M).mat[1, 1]),
            Cnew.size_codomain,
            Cnew.size_domain,
        ) # init new matrix.
        mmIsNonZero = false
        display(mm)
        for (it_ti, ti) in enumerate(eachslice(t, dims = 1))
            for (it_tj, tj) in enumerate(eachslice(t, dims = 1))
                y_test = y - ti + tj
                for (it_yk, yk) in enumerate(eachslice(Ay_old, dims = 1))
                    if isapprox(yk, y_test, rtol = alfa_rtol, atol = alfa_atol)
                        matblock = m_old[it_yk].mat
                        display(matblock)
                        mm[
                           (it_ti-1)*brs+1:it_ti*brs,
                           (it_tj-1)*bcs+1:it_tj*bcs,
                        ] = matblock
                        mmIsNonZero = true
                        break
                    else
                        break
                    end

                end
            end
        end
        if mmIsNonZero
            push!(op, Multiplier(y_new[it_y, :], mm))
        end
    end
    return op
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
        xmax --> xmax+1
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
        for (it_d, sd) in enumerate(eachslice(S.C.Domain, dims = 1))
            for (it_c, sc) in enumerate(eachslice(S.C.Codomain, dims = 1))
                val = real(m.mat[it_d, it_c])
                if abs(val) > threshhold

                    if norm(m.pos) ≈ 0

                    else
                        linecolor := get(
                            ColorSchemes.viridis,
                            -(mmin - val) / (mmax - mmin),
                        ) # linear interpolation
                        label := ""

                        # # draw arrow.
                        p0 = coord_cartesian + sd
                        println("p0")
                        println(p0)
                        p2 = sc
                        println("p2")
                        println(p2)
                        d = -p0 + p2
                        println("d")
                        println(d)

                        mid = p0 .+ d ./ 2
                        println("mid")
                        println(mid)
                        nd = [-d[2], d[1]]
                        p1 = mid + 0.3 * [-d[2], d[1]]#/norm(d)*r
                        println("p1")
                        println(p1)
                        B(t) = (1 - t)^2 * p0 + 2 * (1 - t) * t * p1 + t^2 * p2
                        ## shorten by 10pt

                        #
                        vv = vcat([B(t) for t in range(
                            0,
                            stop = 1,
                            step = 0.01,
                        )]'...)
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
                if abs(val) > threshhold
                    if norm(m.pos) ≈ 0
                        println(val)
                        println(mmin)
                        println(mmax)
                        print(-(mmin - val) / (mmax - mmin))
                        @series begin
                            seriestype --> scatter
                            markercolor := get(
                                ColorSchemes.viridis,
                                -(mmin - val) / (mmax - mmin),
                            ) # linear interpolation
                            #markercolor := get(ColorSchemes.viridis,0.5) # linear interpolation
                            #clim := -5:-4
                            label := ""
                            markershape --> :pentagon
                            markeralpha := 1
                            markersize --> 20 # 10
                            markerstrokewidth --> 0
                            #marker_z = -4:-4
                            pos = coord_cartesian + sd
                            [pos[1]], [pos[2]]
                        end
                    end
                end
            end
        end

    end
end
