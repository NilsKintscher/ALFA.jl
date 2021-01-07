mutable struct CrystalOperator{N,T}
    C::Crystal{N,T} ## dimension L.dim dictates size of Multiplier.pos and dim of structure elements.
    M::SortedSet{Multiplier}# Array{Multiplier,1}
    _CompatibilityCheckOnly::Bool
    function CrystalOperator{N,T}(
        C::Crystal{N,T},
        M::SortedSet{Multiplier},
        _CompatibilityCheckOnly::Bool,
    ) where {N,T<:Union{Float64,Rational}}
        _sanitycheck(C, M)
        new{N,T}(C, M, _CompatibilityCheckOnly)
    end
end


"""
    CrystalOperator(
        C::Crystal{N,T},
        J::UniformScaling,
        _CompatibilityCheckOnly = false,
    ) where {N,T}
    CrystalOperator{N,T}(
        C = nothing,
        M = nothing,
        _CompatibilityCheckOnly = false,
    ) where {N,T<:Union{Float64,Rational}}

Constructs a translationally invariant `C::CrystalOperators{N,T}`
```math
C : \\mathcal{L}(C.L^\\text{Domain}) \\rightarrow \\mathcal{L}(C.L^\\text{Codomain})
```

The actual function definition ``(Cf)(x)`` is saved within `C.M::SortedSet{Multiplier}`:
```math
(Cf)(x) = \\sum_{y \\in C.M} y.\\text{mat} ⋅ f(x+y.\\text{pos}) \\quad \\forall x \\in (C.L.A)\\mathbb{Z}^N
```
# Example
```jldoctest
julia> using ALFA

julia> using LinearAlgebra

julia> ALFA.CrystalOperator{2,Float64}()
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Multiplier: ALFA.Multiplier[]

julia> ALFA.CrystalOperator(ALFA.Crystal{2,Float64}(),3*I)
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Multiplier: 1-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{2}([0, 0], [3])

```
"""
function CrystalOperator(
    C::Crystal{N,T},
    J::UniformScaling,
    _CompatibilityCheckOnly = false,
) where {N,T}
    M = SortedSet{Multiplier}()
    pos = zeros(C.dim)
    mat = Matrix(J, C.size_codomain, C.size_domain)
    push!(M, Multiplier(pos, mat))
    return CrystalOperator{N,T}(C, M, _CompatibilityCheckOnly)
end


function CrystalOperator{N,T}(
    C = nothing,
    M = nothing,
    _CompatibilityCheckOnly = false,
) where {N,T<:Union{Float64,Rational}}
    if C == nothing
        C = Crystal{N,T}()
    end
    if M == nothing
        M = SortedSet{Multiplier}()
    end
    return CrystalOperator{N,T}(C, M, _CompatibilityCheckOnly)
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

"""
    find_multiplier(S::CrystalOperator, pos)

Returns the multiplier at pos if existent (in S.M).

# Example
```jldoctest
julia> using ALFA

julia> S = ALFA.gallery.Laplace()
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Multiplier: 5-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{2}([-1, 0], [1.0])
 ALFA.Multiplier{2}([0, -1], [1.0])
 ALFA.Multiplier{2}([0, 0], [-4.0])
 ALFA.Multiplier{2}([0, 1], [1.0])
 ALFA.Multiplier{2}([1, 0], [1.0])

julia> ALFA.find_multiplier(S, [0, 0])
Position: 2-element StaticArrays.MArray{Tuple{2},Int64,1,2} with indices SOneTo(2):
 0
 0
Multiplier: 1×1 Array{Float64,2}:
 -4.0

```
"""
function find_multiplier(S::CrystalOperator, pos)
    for m in S.M
        if m.pos == pos
            return m
        end
    end
    return nothing
end

"""
    CrystalOperatorCopyWithMultipliers(S::CrystalOperator{N,T}; pos = nothing, idx = nothing) where {N,T}

Returns the CrystalOperator G consisting of the central multiplier of S:
    m_G^{(0)} = p@m_L^{(pos)}@p,
    m_G^{(y)} = 0 if y != pos,

    with P[i,j] = 1 if i=j in idx,
        = 0 else.
"""
function CrystalOperatorCopyWithMultipliers(
    S::CrystalOperator{N,T};
    pos = nothing,
    idx = nothing,
) where {N,T}
    if pos == nothing
        pos = zeros(N)
    end

    m = deepcopy(find_multiplier(S, pos))
    P = zeros(eltype(m.mat), S.C.size_codomain, S.C.size_domain)
    if idx == nothing
        P = I + P
    else
        for j in idx
            P[j, j] = 1
        end
    end

    C = CrystalOperator{N,T}(S.C)
    m.mat =  P * m.mat * P
    push!(C, m)
    return C
end

"""
    CrystalOperatorCopyLowerTriangle(S::CrystalOperator{N,T}; idx = nothing, omega=1.0) where {N,T}

Returns the CrystalOperator G consisting of the central multiplier of S:
    m_G^{(0)} = p@m_L^{(pos)}@p if pos ≦ 0,
    m_G^{(y)} = 0 if pos > 0.

    with P[i,j] = 1 if i=j in idx,
        = 0 else.
"""
function CrystalOperatorCopyLowerTriangle(
    S::CrystalOperator{N,T};
    idx = nothing,
    omega = 1,
    perm = 1:N
) where {N,T}
    pos = zeros(N)

    P = zeros(eltype(first(S.M).mat),S.C.size_codomain, S.C.size_domain)
    if idx == nothing
        P = I + P
    else
        P[idx, idx] .= 1
    end


    C = CrystalOperator{N,T}(S.C)
    for m in S.M
        if m.pos[perm] <= pos[perm]
            mc = deepcopy(m)
            mc.mat = P * mc.mat * P
            if m.pos == pos
                mc.mat = omega * mc.mat
            end

            push!(C, mc)
        end
    end
    return C
end

"""
    Base.push!(S::CrystalOperator, m::Multiplier, add_to_existing = false)

Adds a multiplier to S.M. If there is a multiplier m2 in S.M with m2.pos == m.pos, then
- m2 is replaced if `add_to_existing` == false
- m2.mat+m.mat is the new multiplier at m.pos if `add_to_existing` == true.

# Example
```jldoctest
julia> using ALFA

julia> S = ALFA.gallery.Laplace()
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Multiplier: 5-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{2}([-1, 0], [1.0])
 ALFA.Multiplier{2}([0, -1], [1.0])
 ALFA.Multiplier{2}([0, 0], [-4.0])
 ALFA.Multiplier{2}([0, 1], [1.0])
 ALFA.Multiplier{2}([1, 0], [1.0])

julia> m = ALFA.Multiplier{2}([0, 0], [1])
Position: 2-element StaticArrays.MArray{Tuple{2},Int64,1,2} with indices SOneTo(2):
 0
 0
Multiplier: 1×1 Array{Int64,2}:
 1

julia> push!(S,m)
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Multiplier: 5-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{2}([-1, 0], [1.0])
 ALFA.Multiplier{2}([0, -1], [1.0])
 ALFA.Multiplier{2}([0, 0], [1])
 ALFA.Multiplier{2}([0, 1], [1.0])
 ALFA.Multiplier{2}([1, 0], [1.0])

julia> push!(S,m, true)
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Multiplier: 5-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{2}([-1, 0], [1.0])
 ALFA.Multiplier{2}([0, -1], [1.0])
 ALFA.Multiplier{2}([0, 0], [2])
 ALFA.Multiplier{2}([0, 1], [1.0])
 ALFA.Multiplier{2}([1, 0], [1.0])

```
"""
function Base.push!(S::CrystalOperator, m::Multiplier, add_to_existing = false)
    if add_to_existing
        m_old = find_multiplier(S, m.pos)
        if m_old == nothing
            _sanitycheck(S.C, m)
            push!(S.M, m)
        else
            m_old.mat += m.mat
        end
    else
        _sanitycheck(S.C, m)
        push!(S.M, m)
    end
    return S
end


"""
    normalize(S::CrystalOperator{N,T}) where {N,T}

Normalizes the crystaloperator, i.e., returns an crystaloperator isomorphic to S where all structure elements are shifted in the standard cell ``S.C.L.A\\cdot[0,1)^N`` and sorted lexicographically.

# Example
```jldoctest
julia> using ALFA

julia> S = ALFA.CrystalOperator{1,Float64}(ALFA.Crystal{1,Float64}([1], [-.5], [1.5]))
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [-0.5]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [1.5]
Multiplier: ALFA.Multiplier[]

julia> push!(S, ALFA.Multiplier([0], [-2]))
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [-0.5]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [1.5]
Multiplier: 1-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([0], [-2])

julia> ALFA.normalize(S)
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.5]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.5]
Multiplier: 1-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-2], [-2])

```
"""
function normalize(S::CrystalOperator{N,T}) where {N,T}
    if S.C._IsNormalized
        return S
    end
    (dn, ds, dp) = ShiftIntoStandardCell(S.C.Domain, S.C.L)
    (cn, cs, cp) = ShiftIntoStandardCell(S.C.Codomain, S.C.L)

    return _ModifyStructureElementCore(S, dn, ds, dp, cn, cs, cp)
end

function ChangeStructureElement(S::CrystalOperator{N,T}, Domain, Codomain) where {N,T}
    (dn_from, ds_from, dp_from) = ShiftIntoStandardCell(S.C.Domain, S.C.L)
    (cn_from, cs_from, cp_from) = ShiftIntoStandardCell(S.C.Codomain, S.C.L)

    newC = ALFA.Crystal{N,T}(S.C.L, Domain, Codomain)
    (dn_to, ds_to, dp_to) = ShiftIntoStandardCell(newC.Domain, S.C.L)
    (cn_to, cs_to, cp_to) = ShiftIntoStandardCell(newC.Codomain, S.C.L)

    dn = Domain
    dp = dp_from[invperm(dp_to)]
    ds = ds_from-ds_to[invperm(dp_to)]

    for j in 1:length(Domain)
        @assert dn[j] + S.C.L.A*ds[j] ≈ S.C.Domain[dp[j]] "Domain structure element not compatible?, j= $j"
    end

    cn = Codomain
    cp = cp_from[invperm(cp_to)]
    cs = cs_from-cs_to[invperm(cp_to)]

    for j in 1:length(Codomain)
        @assert cn[j] + S.C.L.A*cs[j] ≈ S.C.Codomain[cp[j]] "Codomain structure element not compatible?, j= $j"
    end

    return _ModifyStructureElementCore(S, dn, ds, dp, cn, cs, cp)
end

function _ModifyStructureElementCore(S::CrystalOperator{N,T}, dn, ds, dp, cn, cs, cp) where {N,T}
    m_old = collect(S.M)
    mattype = typejoin([typeof(mm.mat).parameters[1] for mm in m_old]...)
    if mattype == Any
        mattype = Number
    end
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

    Cnew = Crystal{N,T}(S.C.L.A, dn, cn)
    op = CrystalOperator{N,T}(Cnew)
    op._CompatibilityCheckOnly = S._CompatibilityCheckOnly
    #
    for (y_new, idxset) in SS0
        mat = nothing  # allocate new matrix.
        for (k, i, j) in idxset
            if m_old[k].mat[cp[i], dp[j]] != 0
                if mat == nothing
                    mat = zeros(mattype, size(m_old[1].mat)...) #0 * m_old[1].mat
                end
                mat[i, j] = m_old[k].mat[cp[i], dp[j]]
            end
        end

        if mat != nothing
            push!(op, Multiplier(y_new, mat))
        end
    end
    CleanUp!(op)
    return op
end

"""
    CleanUp!(S::CrystalOperator)

Removes all multipliers m with norm(m.mat) == 0.

# Example
```jldoctest
julia> using ALFA

julia> using LinearAlgebra

julia> S = ALFA.CrystalOperator(ALFA.Crystal{2,Float64}(),0*I)
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Multiplier: 1-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{2}([0, 0], [0])

julia> ALFA.CleanUp!(S)

julia> S
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Multiplier: ALFA.Multiplier[]
```
"""
function CleanUp!(S::CrystalOperator)
    for m in S.M
        if !any(x -> x != 0, m.mat)
            pop!(S.M, m)
        end
    end
end

"""
    wrtLattice(
        S::CrystalOperator{N,T},
        A,
        _CompatibilityCheckOnly = false,
        normalize = true
    ) where {N,T}

Rewrites the crystaloperator S wit hrespect to the translationally invariance A.

# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.gallery.Laplace(N=1)
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [1.0])
 ALFA.Multiplier{1}([0], [-2.0])
 ALFA.Multiplier{1}([1], [1.0])

julia> ALFA.wrtLattice(L, L.C.L.A*2)
Lattice Basis: ALFA.Lattice{1,Float64}([2.0])
Domain: 2-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
 [1.0]
Codomain: 2-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
 [1.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [0.0 1.0; 0.0 0.0])
 ALFA.Multiplier{1}([0], [-2.0 1.0; 1.0 -2.0])
 ALFA.Multiplier{1}([1], [0.0 0.0; 1.0 0.0])
```
"""
function wrtLattice(
    S::CrystalOperator{N,T},
    A,
    _CompatibilityCheckOnly = false,
    normalize = true,
) where {N,T}
    if A isa Lattice
        A = A.A
    end
    t, dH = ElementsInQuotientSpace(
        S.C.A,
        A,
        return_fractional = true,
        return_diag_hnf = true,
    )
    tiMinustj_all = collect(Iterators.product([-x+1:x-1 for x in dH]...))# all possible combinations of t[i]-t[j]

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
    newDomain = [S.C.L.A * x + y for x in t for y in S.C.Domain]

    newCodomain = [S.C.L.A * x + y for x in t for y in S.C.Codomain]

    m_old = collect(S.M)
    mattype = typejoin([typeof(mm.mat).parameters[1] for mm in m_old]...)
    if mattype == Any
        mattype = Number
    end
    y_old = [x.pos for x in S.M]#vcat(transpose([x.pos for x in S.M])...)
    Ay_old = [S.C.L.A * x for x in y_old]#transpose(S.C.L.A * transpose(y_old))


    SS = SortedSet{Array{Int,1}}()
    if T <: Rational
        for j in Iterators.product(y_old, tiMinustj_all)
            y = floor.(inv(A) * (S.C.L.A * (j[1] .+ j[2])))
            # map!(
            #     x -> isapprox(x, round(x), rtol = ALFA_rtol, atol = ALFA_atol) ?
            #         round(x) : floor(x),
            #     y,
            #     y,
            # )
            push!(SS, y)
        end
    else
        for j in Iterators.product(y_old, tiMinustj_all)
            y = A \ (S.C.L.A * (j[1] .+ j[2]))
            map!(
                x ->
                    isapprox(x, round(x), rtol = ALFA_rtol, atol = ALFA_atol) ?
                    round(x) : floor(x),
                y,
                y,
            )
            push!(SS, y)
        end
    end
    y_new = collect(SS)
    Ay_new = [A * x for x in y_new]
    # assign multipliers.

    brs = S.C.size_codomain #block row size
    bcs = S.C.size_domain #block column size

    Cnew = Crystal{N,T}(A, newDomain, newCodomain)
    op = CrystalOperator{N,T}(Cnew)
    op._CompatibilityCheckOnly = _CompatibilityCheckOnly
    for (it_y, y) in enumerate(Ay_new)
        mm = nothing
        for tiMinustj in tiMinustj_all
            y_test = y .- S.C.L.A * [tiMinustj...]
            for (it_yk, yk) in enumerate(Ay_old)
                if isapprox(yk, y_test, rtol = ALFA_rtol, atol = ALFA_atol)
                    matblock = m_old[it_yk].mat
                    if mm == nothing
                        mm =
                            zeros(mattype, Cnew.size_codomain, Cnew.size_domain) # init new matrix.
                    end
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
    if normalize
        op = ALFA.normalize(op)
    end
    CleanUp!(op)
    return op
end

"""
    wrtSameLatticeAndNormalize(A::CrystalOperator, B::CrystalOperator)

Finds a least common multiple translationally invariance C and rewrites both operators A and B wrt C and normalizes the Operators.

# Example
```jldoctest
julia> using ALFA

julia> A = ALFA.gallery.Laplace(N=1)
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [1.0])
 ALFA.Multiplier{1}([0], [-2.0])
 ALFA.Multiplier{1}([1], [1.0])

julia> B = ALFA.gallery.fw_restriction(N=1)
Lattice Basis: ALFA.Lattice{1,Float64}([2.0])
Domain: 2-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
 [1.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 2-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], Rational{Int64}[0//1 1//2])
 ALFA.Multiplier{1}([0], Rational{Int64}[1//1 1//2])

julia> (A2,B2) = ALFA.wrtSameLatticeAndNormalize(A,B);

julia> A2
Lattice Basis: ALFA.Lattice{1,Float64}([2.0])
Domain: 2-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
 [1.0]
Codomain: 2-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
 [1.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [0.0 1.0; 0.0 0.0])
 ALFA.Multiplier{1}([0], [-2.0 1.0; 1.0 -2.0])
 ALFA.Multiplier{1}([1], [0.0 0.0; 1.0 0.0])

julia> B2
Lattice Basis: ALFA.Lattice{1,Float64}([2.0])
Domain: 2-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
 [1.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 2-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], Rational{Int64}[0//1 1//2])
 ALFA.Multiplier{1}([0], Rational{Int64}[1//1 1//2])
```
"""
function wrtSameLatticeAndNormalize(A::CrystalOperator, B::CrystalOperator)
    if A.C.L.A == B.C.L.A
        Anew = A
        Bnew = B
    else
        X = lcm(A.C.L, B.C.L)
        Anew = wrtLattice(A, X)
        Bnew = wrtLattice(B, X)
    end
    Anew = normalize(Anew)
    Bnew = normalize(Bnew)
    return Anew, Bnew
end

function IsApproxEquivalent(A::CrystalOperator, B::CrystalOperator)
    (Anew, Bnew) = wrtSameLatticeAndNormalize(A, B)
    return Anew ≈ Bnew
end



"""
    Base.:/(A::CrystalOperator, b::T) where {T<:Number}

Constructs a translationally invariant `C::CrystalOperators{N,T}`
```math
C/b : \\mathcal{L}(C.L^\\text{Domain}) \\rightarrow \\mathcal{L}(C.L^\\text{Codomain})
```
with
```math
(C/b\\cdot f)(x) = \\frac{1}{b} \\sum_{y \\in C.M} y.\\text{mat} ⋅ f(x+y.\\text{pos}) \\quad \\forall x \\in (C.L.A)\\mathbb{Z}^N
```
# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.gallery.Laplace(N=1)
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [1.0])
 ALFA.Multiplier{1}([0], [-2.0])
 ALFA.Multiplier{1}([1], [1.0])

julia> L/2
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [0.5])
 ALFA.Multiplier{1}([0], [-1.0])
 ALFA.Multiplier{1}([1], [0.5])

```
"""
function Base.:/(A::CrystalOperator, b::T) where {T<:Number}
    return A * (1 / b)
end

"""
    Base.:*(b::T, A::CrystalOperator) where {T<:Number}
    Base.:*(A::CrystalOperator{N,T}, b::S) where {N,T,S<:Number}

Constructs a translationally invariant `C::CrystalOperators{N,T}`
```math
b\\cdot C : \\mathcal{L}(C.L^\\text{Domain}) \\rightarrow \\mathcal{L}(C.L^\\text{Codomain})
```
with
```math
(b\\cdot C\\cdot f)(x) = b\\cdot \\sum_{y \\in C.M} y.\\text{mat} ⋅ f(x+y.\\text{pos}) \\quad \\forall x \\in (C.L.A)\\mathbb{Z}^N
```
# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.gallery.Laplace(N=1)
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [1.0])
 ALFA.Multiplier{1}([0], [-2.0])
 ALFA.Multiplier{1}([1], [1.0])

julia> 2L
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [2.0])
 ALFA.Multiplier{1}([0], [-4.0])
 ALFA.Multiplier{1}([1], [2.0])

```
"""
function Base.:*(b::T, A::CrystalOperator) where {T<:Number}
    return A * b
end
function Base.:*(A::CrystalOperator{N,T}, b::S) where {N,T,S<:Number}
    if A._CompatibilityCheckOnly
        AB = CrystalOperator{N,T}(A.C, nothing, true)
        return AB
    else
        AB = CrystalOperator{N,T}(A.C)
        for am in A.M
            amnew = deepcopy(am)
            amnew.mat *= b
            push!(AB, amnew, true)
        end
        CleanUp!(AB)
        return AB
    end
end

"""
    Base.:*(A::CrystalOperator{N,T}, B::CrystalOperator{N,T}) where {N,T}

If there is a least common multiple translationally invariance of A and B, then
both operators A and B are rewritten wrt this translationally invariance A2 and B2.
After that, it is checked if domain and codomain are compatible. If thats the case,
the crystaloperator A2\\cdot B2 is constructed.
```math
A2\\cdot B2 : \\mathcal{L}(B2.C.L^\\text{Domain}) \\rightarrow \\mathcal{L}(A2.C.L^\\text{Codomain})
```
with
```math
(A2\\cdot B2\\cdot f)(x) = \\sum_{a \\in A2.M, b \\in B2.M} a.\\text{mat}\\cdot b.\\text{mat} ⋅ f(x+a.\\text{pos}+b.\\text{pos}) \\quad \\forall x \\in (A2.C.L.A)\\mathbb{Z}^N.
```
# Example
```jldoctest
julia> using ALFA

julia> A = ALFA.gallery.Laplace(N=1)
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [1.0])
 ALFA.Multiplier{1}([0], [-2.0])
 ALFA.Multiplier{1}([1], [1.0])

julia> R = ALFA.gallery.fw_restriction(N=1)
Lattice Basis: ALFA.Lattice{1,Float64}([2.0])
Domain: 2-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
 [1.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 2-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], Rational{Int64}[0//1 1//2])
 ALFA.Multiplier{1}([0], Rational{Int64}[1//1 1//2])

julia> R*A
Lattice Basis: ALFA.Lattice{1,Float64}([2.0])
Domain: 2-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
 [1.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [0.5 0.0])
 ALFA.Multiplier{1}([0], [-1.0 0.0])
 ALFA.Multiplier{1}([1], [0.5 0.0])

```
"""
function Base.:*(A::CrystalOperator{N,T}, B::CrystalOperator{N,T}) where {N,T}
    if A._CompatibilityCheckOnly || B._CompatibilityCheckOnly
        @assert A.C.L.A == B.C.L.A && isapprox(
            A.C.Domain,
            B.C.Codomain,
            rtol = ALFA_rtol,
            atol = ALFA_atol,
        )
        ABc = Crystal{N,T}(A.C.L, B.C.Domain, A.C.Codomain)
        AB = CrystalOperator{N,T}(ABc, nothing, true)
        return AB
    else
        A, B = wrtSameLatticeAndNormalize(A, B)
        @assert isapprox(
            A.C.Domain,
            B.C.Codomain,
            rtol = ALFA_rtol,
            atol = ALFA_atol,
        )

        ABc = Crystal{N,T}(A.C.L, B.C.Domain, A.C.Codomain)
        AB = CrystalOperator{N,T}(ABc)

        for am in A.M
            for bm in B.M
                y = am.pos + bm.pos
                mat = am.mat * bm.mat
                mult = Multiplier(y, mat)
                push!(AB, mult, true)
            end
        end
        CleanUp!(AB)
        return AB
    end
end


function Base.:-(A::CrystalOperator)
    if A._CompatibilityCheckOnly
        return A
    else
        mA = deepcopy(A)

        for ma in mA.M
            ma.mat *= -1
        end
    end
    return mA
end

function Base.:-(A::CrystalOperator, B::CrystalOperator)
    return A + (-B)
end

"""
    Base.:+(A::CrystalOperator{N,T}, B::CrystalOperator{N,T}) where {N,T}

If there is a least common multiple translationally invariance of A and B, then
both operators A and B are rewritten wrt this translationally invariance A2 and B2.
After that, it is checked if domain and codomain are compatible. If thats the case,
the crystaloperator A2+t B2 is constructed.
```math
A2+B2 : \\mathcal{L}(A2.C.L^\\text{Domain}) \\rightarrow \\mathcal{L}(A2.C.L^\\text{Codomain})
```
with
```math
(A2+B2\\cdot f)(x) = (A2\\cdot f)(x) + (B2\\cdot f)(x) \\quad \\forall x \\in (A2.C.L.A)\\mathbb{Z}^N.
```
# Example
```jldoctest
julia> using ALFA

julia> A = ALFA.gallery.Laplace(N=1)
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [1.0])
 ALFA.Multiplier{1}([0], [-2.0])
 ALFA.Multiplier{1}([1], [1.0])

julia> A+A
Lattice Basis: ALFA.Lattice{1,Float64}([1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]
Multiplier: 3-element Array{ALFA.Multiplier,1}:
 ALFA.Multiplier{1}([-1], [2.0])
 ALFA.Multiplier{1}([0], [-4.0])
 ALFA.Multiplier{1}([1], [2.0])

```
"""
function Base.:+(A::CrystalOperator{N,T}, B::CrystalOperator{N,T}) where {N,T}
    if A._CompatibilityCheckOnly || B._CompatibilityCheckOnly
        @assert A.C.L.A == B.C.L.A &&
                isapprox(
                    A.C.Domain,
                    B.C.Domain,
                    rtol = ALFA_rtol,
                    atol = ALFA_atol,
                ) &&
                isapprox(
                    A.C.Codomain,
                    B.C.Codomain,
                    rtol = ALFA_rtol,
                    atol = ALFA_atol,
                )
        AB = CrystalOperator{N,T}(A.C, nothing, true)
        return AB
    else
        A, B = wrtSameLatticeAndNormalize(A, B)
        @assert isapprox(
            A.C.Domain,
            B.C.Domain,
            rtol = ALFA_rtol,
            atol = ALFA_atol,
        ) && isapprox(
            A.C.Codomain,
            B.C.Codomain,
            rtol = ALFA_rtol,
            atol = ALFA_atol,
        )

        AB = deepcopy(A)


        for bm in B.M
            push!(AB, bm, true)
        end

        CleanUp!(AB)
        return AB
    end
end


function Base.transpose(A::CrystalOperator{N,T}) where {N,T}
    if A._CompatibilityCheckOnly
        tA = CrystalOperator{N,T}(
            Crystal{N,T}(A.C.L, A.C.Codomain, A.C.Domain),
            nothing,
            true,
        )
        return tA
    else
        tA = CrystalOperator{N,T}(Crystal{N,T}(A.C.L, A.C.Codomain, A.C.Domain))

        for ma in A.M
            m = Multiplier(-ma.pos, transpose(ma.mat))
            push!(tA, m)
        end
    end
    return tA
end

function Base.adjoint(A::CrystalOperator{N,T}) where {N,T}
    if A._CompatibilityCheckOnly
        tA = CrystalOperator{N,T}(
            Crystal{N,T}(A.C.L, A.C.Codomain, A.C.Domain),
            nothing,
            true,
        )
        return tA
    else
        tA = CrystalOperator{N,T}(Crystal{N,T}(A.C.L, A.C.Codomain, A.C.Domain))

        for ma in A.M
            m = Multiplier(-ma.pos, adjoint(ma.mat))
            push!(tA, m)
        end
    end
    return tA
end

function LinearAlgebra.pinv(A::CrystalOperator)
    @assert A._CompatibilityCheckOnly
    return transpose(A)
end
function LinearAlgebra.inv(A::CrystalOperator)
    @assert A._CompatibilityCheckOnly
    return transpose(A)
end

function Base.:^(A::CrystalOperator, p::Int)
    @assert p > 0
    Ac = normalize(A)
    return prod(Ac for _ = 1:p)
end

function Base.:+(J::UniformScaling, A::CrystalOperator)
    return A + J
end
function Base.:+(A::CrystalOperator, J::UniformScaling)
    Ac = deepcopy(A)
    pos = zeros(A.C.dim)
    mat = Matrix(J, A.C.size_codomain, A.C.size_domain)
    push!(Ac, Multiplier(pos, mat), true)
    CleanUp!(Ac)
    return Ac
end

function Base.:-(A::CrystalOperator, J::UniformScaling)
    return A + (-J)
end
function Base.:-(J::UniformScaling, A::CrystalOperator)
    return J + (-A)
end



function Base.:(==)(A::CrystalOperator, B::CrystalOperator)
    if A.C == B.C && length(A.M) == length(B.M)
        for (ma, mb) in zip(A.M, B.M)
            if ma != mb
                return false
            end
        end
        return true
    else
        return false
    end
end



function Base.:(≈)(A::CrystalOperator, B::CrystalOperator)
    if A.C ≈ B.C && length(A.M) == length(B.M)
        for (ma, mb) in zip(A.M, B.M)
            if !(ma ≈ mb)
                return false
            end
        end
        return true
    else
        return false
    end
end


"""
    construct_matrix(x::ALFA.CrystalOperator,wrt::ALFA.Lattice)

Constructs a matrix from a CrystalOperator with respect to a given Lattice.
"""
function construct_matrix(x::ALFA.CrystalOperator,wrt::ALFA.Lattice)
    xc = ALFA.wrtLattice(x, wrt)
    xc = ALFA.normalize(xc)
    xmat = sum(x.mat for x in xc.M);
    return xmat
end
