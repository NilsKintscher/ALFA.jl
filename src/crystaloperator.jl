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

function find_multiplier(S::CrystalOperator, pos)
    for m in S.M
        if m.pos == pos
            return m
        end
    end
    return nothing
end

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



function normalize(S::CrystalOperator{N,T}) where {N,T}
    if S.C._IsNormalized
        return S
    end
    (dn, ds, dp) = ShiftIntoUnitCell(S.C.Domain, S.C.L)
    (cn, cs, cp) = ShiftIntoUnitCell(S.C.Codomain, S.C.L)


    m_old = collect(S.M)
    mattype = typejoin([typeof(mm.mat).parameters[1] for mm in m_old]...)
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

function CleanUp!(S::CrystalOperator)
    for m in S.M
        if !any(x -> x != 0, m.mat)
            pop!(S.M, m)
        end
    end
end

function wrtLattice(
    S::CrystalOperator{N,T},
    A,
    _CompatibilityCheckOnly = false,
) where {N,T} ### A::Matrix
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
    y_old = [x.pos for x in S.M]#vcat(transpose([x.pos for x in S.M])...)
    Ay_old = [S.C.L.A * x for x in y_old]#transpose(S.C.L.A * transpose(y_old))


    SS = SortedSet{Array{Int,1}}()
    if T <: Rational
        for j in Iterators.product(y_old, tiMinustj_all)
            y = floor.(inv(A) * (S.C.L.A * (j[1] .+ j[2])))
            # map!(
            #     x -> isapprox(x, round(x), rtol = alfa_rtol, atol = alfa_atol) ?
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
                    isapprox(x, round(x), rtol = alfa_rtol, atol = alfa_atol) ?
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
                if isapprox(yk, y_test, rtol = alfa_rtol, atol = alfa_atol)
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
    CleanUp!(op)
    return op
end

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



function Base.:/(A::CrystalOperator, b::T) where {T<:Number}
    return A * (1 / b)
end
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

function Base.:*(A::CrystalOperator{N,T}, B::CrystalOperator{N,T}) where {N,T}
    if A._CompatibilityCheckOnly || B._CompatibilityCheckOnly
        @assert A.C.L.A == B.C.L.A && isapprox(
            A.C.Domain,
            B.C.Codomain,
            rtol = alfa_rtol,
            atol = alfa_atol,
        )
        ABc = Crystal{N,T}(A.C.L, B.C.Domain, A.C.Codomain)
        AB = CrystalOperator{N,T}(ABc, nothing, true)
        return AB
    else
        A, B = wrtSameLatticeAndNormalize(A, B)
        @assert isapprox(
            A.C.Domain,
            B.C.Codomain,
            rtol = alfa_rtol,
            atol = alfa_atol,
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

function Base.:+(A::CrystalOperator{N,T}, B::CrystalOperator{N,T}) where {N,T}
    if A._CompatibilityCheckOnly || B._CompatibilityCheckOnly
        @assert A.C.L.A == B.C.L.A &&
                isapprox(
                    A.C.Domain,
                    B.C.Domain,
                    rtol = alfa_rtol,
                    atol = alfa_atol,
                ) &&
                isapprox(
                    A.C.Codomain,
                    B.C.Codomain,
                    rtol = alfa_rtol,
                    atol = alfa_atol,
                )
        AB = CrystalOperator{N,T}(A.C, nothing, true)
        return AB
    else
        A, B = wrtSameLatticeAndNormalize(A, B)
        @assert isapprox(
            A.C.Domain,
            B.C.Domain,
            rtol = alfa_rtol,
            atol = alfa_atol,
        ) && isapprox(
            A.C.Codomain,
            B.C.Codomain,
            rtol = alfa_rtol,
            atol = alfa_atol,
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
