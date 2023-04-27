

struct Lattice{N,T}
    A::MMatrix{N,N,T}
    function Lattice{N,T}(
        A::MMatrix{N,N,T},
    ) where {N,T<:Union{Float64,Rational}}
        @assert !isapprox(det(A), 0) "Basis must be nonsingular"
        new{N,T}(A)
    end
end
# \{∑_{j=1}^n a_j z_j \ : \ z_j ∈ ℤ\}.
"""
    Lattice{N,T}(A=nothing)

Construct a `Lattice{N,T}` with basis ``A``, i.e., it represents the set
```math
Aℤ^N = \\left\\{ ∑_{i=0}^n z_j ⋅ a_j  : z_j ∈ ℤ \\right\\},
```
where ``a_j`` denotes the ``j``th column of ``A``. The matrix ``A`` must be square and nonsingular:
- ``A`` is called the *lattice basis*.
- ``a_j`` are the *primite vectors* of the lattice.
`T<:Union{Float64,Rational}` represents the datatype of the entries of the basis ``a_{ij}``.

In case of A==nothing, the identity I of size N is used.

See Definition 2.1 in [1].

# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.Lattice{3,Float64}()
ALFA.Lattice{3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])

julia> L = ALFA.Lattice{2,Rational{BigInt}}()
ALFA.Lattice{2, Rational{BigInt}}(Rational{BigInt}[1//1 0//1; 0//1 1//1])

julia> L = ALFA.Lattice{2,Rational{BigInt}}([1 2; 3 4])
ALFA.Lattice{2, Rational{BigInt}}(Rational{BigInt}[1//1 2//1; 3//1 4//1])

julia> L = ALFA.Lattice{2,Rational{BigInt}}([1 2; 1 2])
ERROR: AssertionError: Basis must be nonsingular
```
"""
function Lattice{N,T}(A = nothing) where {N,T<:Union{Float64,Rational}}
    if A == nothing
        A = MMatrix{N,N,T}(I) # identity matrix.
    elseif A isa Matrix
        @assert size(A, 1) == size(A, 2) "Matrix must be square."
        A = MMatrix{N,N,T}(A)
    elseif A isa Real || A isa UniformScaling
        A = MMatrix{N,N,T}(A)
    elseif A isa Vector
        A = MMatrix{N,N,T}(A...)
    else
        A = convert(MMatrix{N,N,T}, A)
    end
    Lattice{N,T}(A)
end
"""
    Base.size(L::Lattice)

Returns the dimension of the lattice basis: size(L.A,1).

# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.Lattice{2, Float64}([1 2; 3 4])
ALFA.Lattice{2, Float64}([1.0 2.0; 3.0 4.0])
julia> size(L) == size(L.A) == (2,2)
true
```
"""
Base.size(L::Lattice) = size(L.A)

"""
    Base.getindex(L::Lattice, y...)

simply wrapping getindex(L.A, y...).

# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.Lattice{2, Float64}([1 2; 3 4])
ALFA.Lattice{2, Float64}([1.0 2.0; 3.0 4.0])
julia> L[1,2] == L.A[1,2]
true
```
"""
Base.getindex(L::Lattice, y...) = getindex(L.A, y...)

"""
    Base.getproperty(L::Lattice, sym::Symbol)

Get properties of Lattice{N,T}. Let ``A=``L.A, i.e., ``A ∈ T^{N × N}``, then
- L.n and L.dim return the dimension N.
- L.iA returns ``A^{-1}``, i.e., the *inverse* of the lattice basis
- L.dA returns ``A^{-T}``, i.e., the basis of the *dual* lattice.

# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.Lattice{2, Float64}([1 2; 3 4])
ALFA.Lattice{2, Float64}([1.0 2.0; 3.0 4.0])

julia> L.n == L.dim == typeof(L).parameters[1] == 2
true
julia> L[2,1]
3.0
julia> L.iA == inv(L.A)
true
julia> L.dA == transpose(inv(L.A))
true
```
"""
function Base.getproperty(L::Lattice, sym::Symbol)
    if sym == :n || sym == :dim
        typeof(L).parameters[1]
    elseif sym == :iA
        inv(L.A)
    elseif sym == :dA
        transpose(inv(L.A))
    else
        # fallback to getfield
        getfield(L, sym)
    end
end


"""
    Base.lcm(A::MArray{X,T},B::MArray{X,T}...)
    Base.lcm(A::MArray{X,T}, B::MArray{X,T}...) where {X,T<:Rational}

    Base.lcm(A::Lattice{N,T}, B::Lattice{N,T}...) where {N,T}
    Base.lcm(A::MArray{X,T}, B::MArray{X,T}...) where {X,T<:Rational}



Returns the least common multiple of ``A`` and ``B`` (or more), i.e. a sub-lattice ``C``, that is ``C ⊂ A`` and ``C ⊂ B`` with ``|\\det(C)|`` as small as possible.

See Algorithm B.7, Theorem 5.1 in [1].
# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.Lattice{2, Float64}([1 2; 3 4])
ALFA.Lattice{2, Float64}([1.0 2.0; 3.0 4.0])

julia> ALFA.lcm(L) == L
true

julia> ALFA.lcm(L.A) == L.A
true

julia> ALFA.lcm(ALFA.Lattice{1,Float64}([2]), ALFA.Lattice{1,Float64}([3]), ALFA.Lattice{1,Float64}([6]))
ALFA.Lattice{1, Float64}([6.0;;])

julia> ALFA.lcm(ALFA.Lattice{2,Rational{BigInt}}([1 1; -1 1]), ALFA.Lattice{2,Rational{BigInt}}([1 2; 2 1]))
ALFA.Lattice{2, Rational{BigInt}}(Rational{BigInt}[-3//1 1//1; -3//1 -1//1])
```
"""
function Base.lcm(A::MArray{X,T}, B::MArray{X,T}) where {X,T<:Real}
    M0 = inv(A) * B
    for r in range(1, stop = 50)
        rM = r * M0
        rMr = round.(rM)
        if isapprox(rM, rMr, rtol = ALFA_rtol, atol = ALFA_atol)
            # compute SNF
            S, U, V = snf_with_transform(rMr)
            N = diagm(r ./ gcd.(diag(S), r))
            C = B * ALFA.lll(V * N)
            return convert(typeof(A), C)
            break
        end
    end
end

function Base.lcm(A::MArray{X,T}, B::MArray{X,T}) where {X,T<:Rational}
    M0 = inv(A) * B
    r = lcm(denominator.(M0)...)
    rM = r * M0
    S, U, V = snf_with_transform(rM)
    N = diagm(r ./ gcd.(diag(S), r))
    C = B * ALFA.lll(V * N)
    return convert(typeof(A), C)
end

function Base.lcm(A::Lattice{N,T}, B::Lattice{N,T}) where {N,T}
    C = lcm(A.A, B.A)
    return Lattice{N,T}(C)
end

Base.lcm(A::Lattice{N,T}, B::Lattice{N,T}...) where {N,T} =
    Base.lcm(A, Base.lcm(B...))
Base.lcm(A::MArray{X,T}, B::MArray{X,T}...) where {X,T} =
    Base.lcm(A, Base.lcm(B...))
Base.lcm(A::MArray{X,T}) where {X,T} = A
Base.lcm(A::Lattice{N,T}) where {N,T} = A

"""
    ElementsInQuotientSpace(
        A::Union{Matrix{T},MMatrix{N,N,T}},
        B::Union{Matrix{T},MMatrix{N,N,T}};
        return_diag_hnf::Bool = false,
        return_fractional::Bool = false,
    ) where {N,T}



Returns all lattice points of the lattice generated by ``A`` found in the primitive cell of B, i.e.,
```math
T_{A,B}=\\{x : x ∈ A \\mathbb{Z}^N \\cap B[0,1)^N \\}
```

See Algorithm B.3, Theorem 2.3 in [1].
# Example
```jldoctest
julia> using ALFA

julia> A = ALFA.Lattice{2, Float64}([1 0; 0 1]);

julia> B = ALFA.Lattice{2, Float64}([2 0; 0 2]);

julia> ALFA.ElementsInQuotientSpace(A.A,B.A)
4-element Vector{StaticArraysCore.SVector{2, Float64}}:
 [0.0, 0.0]
 [1.0, 0.0]
 [0.0, 1.0]
 [1.0, 1.0]
```
"""
function ElementsInQuotientSpace(
    A::Union{Matrix{T},MMatrix{N,N,T}},
    B::Union{Matrix{T},MMatrix{N,N,T}};
    return_diag_hnf::Bool = false,
    return_fractional::Bool = false,
) where {N,T}

    M = inv(A) * B
    Mr = round.(M)
    @assert isapprox(M, Mr, rtol = ALFA_rtol, atol = ALFA_atol) "A must be a sublattice of B, i.e. A\\B must be integral."
    H = convert(Matrix{Int}, hnf(Mr))
    dH = diag(H)
    m = prod(dH)
    J = Iterators.product([0:x-1 for x in dH]...)

    if return_fractional
        if T == BigInt
            s = [SVector{length(dH),BigInt}(x) for x in J]
        else
            s = [SVector{length(dH),Int}(x) for x in J]
        end
    else
        s = [SVector{length(dH),T}(A * [x...]) for x in J]
    end
    if return_diag_hnf
        return vec(s), dH
    else
        return vec(s)
    end
end

"""
    ShiftIntoStandardCell(s, A::Union{Matrix{T},MMatrix{N,N,T}}) where {N,T}
    ShiftIntoStandardCell(s, A::Union{Matrix{T},MMatrix{N,N,T}}) where {N,T<:Rational}
    ShiftIntoStandardCell(s, A::Lattice)

Shifts all elements s[i] into the primitive cell ``A[0,1)^N`` and sort the entries lexicographically.
The function returns t, y, p, such that
```math
t_j + A y_j = s_{p(j)}
```
- t is the shifted vector s, i.e., ``A^{-1}t_j ∈ [0,1)^N`` for all j.
- p is the permutation with respect to the input s.
- y corresponds to the shift in fractional coordinates.

See Definition 2.2, 5.2  in [1].
# Example
```jldoctest
julia> using ALFA

julia> using LinearAlgebra

julia> L = ALFA.Lattice{2, Float64}([1 0; 0 1]);

julia> s = [[-1/2,0], [1/4,2]]
2-element Vector{Vector{Float64}}:
 [-0.5, 0.0]
 [0.25, 2.0]

julia> (t, y, p) = ALFA.ShiftIntoStandardCell(s,L);

julia> t
2-element Vector{StaticArraysCore.MVector{2, Float64}}:
 [0.25, 0.0]
 [0.5, 0.0]

julia> y
2-element Vector{StaticArraysCore.MVector{2, Float64}}:
 [0.0, 2.0]
 [-1.0, 0.0]

julia> p
2-element Vector{Int64}:
 2
 1

julia> [t[j] + L.A*y[j] - s[p[j]] for j in [1,2]]
2-element Vector{StaticArraysCore.MVector{2, Float64}}:
 [0.0, 0.0]
 [0.0, 0.0]
```
"""
function ShiftIntoStandardCell(s, A::Union{Matrix{T},MMatrix{N,N,T}}) where {N,T}
    y = [MVector{N}(A \ x) for x in s]
    #@show y
    for yy in y
        map!(
            x -> isapprox(x, round(x), rtol = ALFA_rtol, atol = ALFA_atol) ?
                round(x) : floor(x),
            yy,
            yy,
        ) # remove non-integral part ∈(0,1)
    end
    #y = [floor.(x) for x in y]
    Ay = [A * x for x in y]
    #transpose(A * transpose(y))
    s = s - Ay
    # sort lexicographically
    p = sortperm(s, alg = Base.Sort.DEFAULT_STABLE, lt = islessapprox)
    s = s[p]
    y = y[p]
    return s, y, p
end

function ShiftIntoStandardCell(
    s,
    A::Union{Matrix{T},MMatrix{N,N,T}},
) where {N,T<:Rational} # s vector of SVector

    y = [inv(A) * x for x in s]
    y = [floor.(x) for x in y]
    Ay = [A * x for x in y]
    s = s - Ay

    # sort lexicographically
    p = sortperm(s, alg = Base.Sort.DEFAULT_STABLE, lt = isless)
    s = s[p]
    y = y[p]
    return s, y, p # We now have s in the primitive cell of A (A*[0,1)^dim) and s is lexicographically ordered
end

function ShiftIntoStandardCell(s, A::Lattice)
    return ShiftIntoStandardCell(s, A.A)
end

"""
    CheckIfNormal(s, A)

Checks if s is found within the primtive cell of A, i.e., ``A^{-1}s_j ∈ [0,1)^N`` and if s is sorted lexicographically.

See Definition 5.4 in [1].
# Example
```jldoctest
julia> using ALFA

julia> using LinearAlgebra

julia> L = ALFA.Lattice{2, Float64}([1 0; 0 1]);

julia> s = [[-1/2,0], [1/4,2]]
2-element Vector{Vector{Float64}}:
 [-0.5, 0.0]
 [0.25, 2.0]

julia>  ALFA.CheckIfNormal(s,L)
false

julia> t = [[1/4,0], [1/2,0]]
2-element Vector{Vector{Float64}}:
 [0.25, 0.0]
 [0.5, 0.0]

julia> ALFA.CheckIfNormal(t,L)
true
```
"""
function CheckIfNormal(s, A)
    (n, s, p) = ShiftIntoStandardCell(s, A)

    @assert n == unique(n) "The points of the structure element must be unique."
    if s == 0 * s && p == sort(p)
        return true
    end
    return false
end


function Base.:(==)(A::Lattice, B::Lattice)
    return A.A == B.A
end


function Base.:(≈)(A::Lattice, B::Lattice)
    return A.A ≈ B.A
end
