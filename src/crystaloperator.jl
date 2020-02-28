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

function CrystalOperator(C = nothing, M = nothing, _CompatibilityCheckOnly = false)
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
    (dn,ds,dp) = ShiftIntoUnitCell(S.C.Domain, S.C.L)
    (cn,cs,cp) = ShiftIntoUnitCell(S.C.Codomain, S.C.L)

    # ...

end

function wrtLattice(S::CrystalOperator, A::Matrix)
    #### TODO: Test this function, when Plot function exists. Need to construct test cases with alfa.py.
    t = ElementsInQuotientSpace(S.C.A, A, fractional=true)
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
    Ay_old = transpose(S.C.L.A*transpose(y_old))
    # find all unique combinations of floor(A\S.C.L.A*(y_old[i] + t[j] - t[k]))
    SS = SortedSet{Array{Int,1}}()
    for j in Iterators.product(eachslice(y_old,dims=1), eachslice(t,dims=1), eachslice(t,dims=1))
        y = A\(S.C.L.A*(j[1]+j[2]-j[3]))
        map!(x -> isapprox(x,round(x), rtol=alfa_rtol, atol=alfa_atol) ? round(x) : floor(x), y, y)
        push!(SS, y)
    end
    y_new = vcat(transpose(collect(SS))...)

    Ay_new = transpose(A*transpose(y_new)) # convert to cartesian coordinate of original lattice.

    # get coordinate of new lattice
    print(Ay_new)
    # assign multipliers.

    brs = S.C.size_codomain #block row size
    bcs = S.C.size_domain #block column size

    Cnew = Crystal(A, newDomain, newCodomain)
    op = CrystalOperator(Cnew)
    for (it_y, y) in enumerate(eachslice(Ay_new,dims=1))
        mm = zeros(typeof(first(S.M).mat[1,1]), Cnew.size_codomain, Cnew.size_domain) # init new matrix.
        mmIsNonZero = false
        display(mm)
        for (it_ti,ti) in enumerate(eachslice(t,dims=1))
            for (it_tj,tj) in enumerate(eachslice(t,dims=1))
                y_test = y - ti + tj
                for (it_yk, yk) in enumerate(eachslice(Ay_old,dims=1))
                    if isapprox(yk, y_test, rtol=alfa_rtol, atol= alfa_atol)
                        matblock = m_old[it_yk].mat
                        display(matblock)
                        mm[(it_ti-1)*brs+1:it_ti*brs, (it_tj-1)*bcs+1:it_tj*bcs] = matblock
                        mmIsNonZero = true
                        break
                    else
                        break
                    end

                end
            end
        end
        if mmIsNonZero
            push!(op, Multiplier(y_new[it_y,:], mm))
        end
    end
    return op
end
