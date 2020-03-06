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

function symbol(S::CrystalOperator, k)
    if length(S.M) == 0
        mat = zeros(Int, S.C.size_codomain, S.C.size_domain)
    else
        mat = zeros(eltype(first(S.M).mat), S.C.size_codomain, S.C.size_domain)
        for m in S.M
            mat += m.mat * exp(im * 2π * dot((S.C.A * m.pos), k))
        end
    end
    return mat
end

function eigvals(S::CrystalOperator, k; by=abs)
    symb = symbol(S, k)
    ev = LinearAlgebra.eigvals(symb)
    ev = sort(ev, by=by)
    return ev
end

function eigen(S::CrystalOperator, k; by=abs)
    symb = symbol(S, k)
    ev = LinearAlgebra.eigen(symb)
    p = sortperm(ev.values, by=by)
    return LinearAlgebra.Eigen(ev.values[p], ev.vectors[:,p])
end

function normalize(S::CrystalOperator)

    (dn, ds, dp) = ShiftIntoUnitCell(S.C.Domain, S.C.L)
    (cn, cs, cp) = ShiftIntoUnitCell(S.C.Codomain, S.C.L)


    m_old = collect(S.M)
    y_old = [x.pos for x in S.M]#vcat(transpose([x.pos for x in S.M])...)
    # find all combinations of the shifts -ds+cs
    SS0 = SortedDict{Array{Int,1},Array{Tuple{Int,Int,Int},1}}() # find more efficient way.
    for (k, yk) in enumerate(y_old) #eachslice(y_old, dims = 1))
        for (i, csi) in enumerate(cs) # eachslice(cs, dims = 1))
            for (j, dsj) in enumerate(ds) #eachslice(ds, dims = 1))
                push!(
                    get!(SS0, yk + dsj - csi, Array{Tuple{Int64,Int64},1}()),
                    (k, i, j),
                )
            end
        end
    end

    Cnew = Crystal(S.C.L.A, dn, cn)
    op = CrystalOperator(Cnew)
    #
    for (y_new, idxset) in SS0
        mat = nothing  # allocate new matrix.
        for (k, i, j) in idxset
            if m_old[k].mat[cp[i], cp[j]] != 0
                if mat == nothing
                    mat = 0 * m_old[1].mat
                end
                mat[i, j] = m_old[k].mat[cp[i], cp[j]]
            end
        end

        if mat != nothing
            push!(op, Multiplier(y_new, mat))
        end
    end
    return op
end


function wrtLattice(S::CrystalOperator, A) ### A::Matrix
    #### TODO: Test this function, when Plot function exists. Need to construct test cases with alfa.py.

    t, dH = ElementsInQuotientSpace(
        S.C.A,
        A,
        return_fractional = true,
        return_diag_hnf = true,
    )
    tiMinustj_all = collect(Iterators.product([-x+1:x-1 for x in dH]...)) # all possible combinations of t[i]-t[j]

    # helper functions for t and index calculation of t
    function RangeOfTi(tiMinustj)
        # get all possible combinations of ti, such that tiMinustj = ti - tj with ti and tj in t
        return Iterators.product([
            max(0, kh[1]):min((kh[2] - 1), (kh[2] - 1) + kh[1])
            for kh in zip(tiMinustj, dH)
        ]...)
    end
    dHprod = [prod(dH[1:i-1]) for i in eachindex(dH)]
    function lookup_idxij(ti, tiMinustj)
        # given t[i] and t[i]-t[j], returns the tuple (i,j).
        tj = ti .- tiMinustj
        return sum(x * y for (x, y) in zip(ti, dHprod)) + 1,
        sum(x * y for (x, y) in zip(tj, dHprod)) + 1
    end
    #
    newDomain = [
        S.C.L.A*x + y for x in t
        for y in S.C.Domain
    ]

    newCodomain = [
        S.C.L.A*x + y for x in t
        for y in S.C.Codomain
    ]

    m_old = collect(S.M)
    y_old = [x.pos for x in S.M]#vcat(transpose([x.pos for x in S.M])...)
    Ay_old = [S.C.L.A * x for x in y_old]#transpose(S.C.L.A * transpose(y_old))


    SS = SortedSet{Array{Int,1}}()
    for j in Iterators.product(y_old, tiMinustj_all)
        y = A \ (S.C.L.A * (j[1] .+ j[2]))
        map!(
            x -> isapprox(x, round(x), rtol = alfa_rtol, atol = alfa_atol) ?
                round(x) : floor(x),
            y,
            y,
        )
        push!(SS, y)
    end
    y_new = collect(SS) #vcat(transpose(collect(SS))...)
    #println("y_new", y_new)
    Ay_new = [A*x for x in y_new] #transpose(A * transpose(y_new)) # convert to cartesian coordinate of original lattice.
    #Ay_new = y_new
    # get coordinate of new lattice
    #println("Ay_new: ", Ay_new)
    #println("Ay_old:", Ay_old)
    # assign multipliers.

    brs = S.C.size_codomain #block row size
    bcs = S.C.size_domain #block column size

    Cnew = Crystal(A, newDomain, newCodomain)
    op = CrystalOperator(Cnew)
    for (it_y, y) in enumerate(Ay_new) # eachslice(Ay_new, dims = 1))
        #println("y_new: ",y)
        mm = nothing
        for tiMinustj in tiMinustj_all
            #for (tdiff, ss_ij) in SS0

            #for (it_ti, ti) in enumerate(eachslice(t, dims = 1))
            #for (it_tj, tj) in enumerate(eachslice(t, dims = 1))
            #y_test = y - ti + tj
            #y_test = y - tdiff
            y_test = y .- S.C.L.A*[tiMinustj...]
            #println("y_test  $y_test , it (i,j) = ($it_ti, $it_tj)")
            for (it_yk, yk) in enumerate(Ay_old)
                #print("yk $yk")
                if isapprox(yk, y_test, rtol = alfa_rtol, atol = alfa_atol)
                    #println("#####################TRUE : yk:   $yk")
                    matblock = m_old[it_yk].mat
                    #println("matblock: $matblock")
                    if mm == nothing
                        mm = zeros(
                            eltype(first(S.M).mat),
                            Cnew.size_codomain,
                            Cnew.size_domain,
                        ) # init new matrix.
                    end
                    #for (it_ti, it_tj) in ss_ij
                    for (it_ti, it_tj) in [
                        lookup_idxij(x, tiMinustj)
                        for x in RangeOfTi(tiMinustj)
                    ]
                        mm[
                            (it_ti-1)*brs+1:it_ti*brs,
                            (it_tj-1)*bcs+1:it_tj*bcs,
                        ] = matblock
                    end
                    break
                end

                #end
            end
        end
        if mm != nothing
            push!(op, Multiplier(y_new[it_y], mm))
        end
    end
    return op
end
