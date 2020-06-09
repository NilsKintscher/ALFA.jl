
mutable struct Crystal{N, T}
    L::Lattice{N, T}
    Domain::Vector{SVector{N,T}}   #Matrix # m × dim
    Codomain::Vector{SVector{N,T}} # m × dim
    _IsNormalized::Bool
    function Crystal{N, T}(
        L::Lattice{N, T},
        Domain::Vector{SVector{N,T}},
        Codomain::Vector{SVector{N,T}},
        _IsNormalized::Bool
    ) where {N, T<:Number}#
        new{N, T}(L, Domain, Codomain, _IsNormalized)
    end
end

"""
    Crystal{N,T}(L = nothing, Domain = nothing, Codomain = nothing) where {N,T<:Union{Float64, Rational}}

Constructs a `Crystal{N,T}` which respresents the domain and codomain of a CrystalOperator.
It consists of a L::Lattice{N,T} and the structure elements Domain::Vector{SVector{N,T}} and Codomain::Vector{SVector{N,T}}.

This struct describes the set
```math
Aℤ^N+s = \\left\\{ ∑_{i=0}^n z_j ⋅ a_j + (s_1,s_2,\\ldots,s_m)  : z_j ∈ ℤ  \\right\\},
```
where ``s=(s_1,s_2,\\ldots,s_m) ∈ \\{\\text{Domain,Codomain}\\}`` is the *structure element*, ``A`` the lattice basis L.A.

- In case of L==nothing, the identity I of size N is used.
- In case of Domain==nothing, the origin `` 0 \\in \\mathbb{R}^N`` is used.
- In case Codomain==nothing, Domain is used.

# Example
```jldoctest
julia> using ALFA

julia> ALFA.Crystal{2,Float64}()
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
```
"""
function Crystal{N,T}(L = nothing, Domain = nothing, Codomain = nothing) where {N,T<:Union{Float64, Rational}}
    if L == nothing
        L = Lattice{N,T}()
    elseif !isa(L, Lattice) #typeof(L) != Lattice
        L = Lattice{N,T}(L)
    end

    #N = L.dim# typeof(L).parameters[1] # dimension

    if Domain == nothing # put a point at the origin.
        Domain = [zeros(SVector{N,T})]
    elseif typeof(Domain) <: Vector{Vector{T}} where {T<:Real} ||
           Domain isa Vector{Matrix{T}} where {T<:Real} #  turn Vector of Vector into Vector of SVector
        Domain = [SVector{N,T}(x) for x in Domain]
    elseif Domain isa Vector{<:Number} # turn vector into vector of SVector.
        if N == 1
            Domain = [SVector{N,T}(x) for x in Domain]
        else
            Domain = [SVector{N,T}(Domain)]
        end
    elseif Domain isa Matrix{<:Number} # turn a matrix rowwise into a vector of SVector.
        Domain = [SVector{N,T}(x) for x in eachrow(Domain)]
    end

    if Codomain == nothing # put a point at the origin.
        Codomain = Domain
    elseif typeof(Codomain) <: Vector{Vector{T}} where {T<:Real} ||
           Codomain isa Vector{Matrix{T}} where {T<:Real}# turn Vector into Matrix with 1 Column
        Codomain = [SVector{N,T}(x) for x in Codomain]
    elseif Codomain isa Vector{<:Number} # turn vector into vector of SVector.
        if N == 1
            Codomain = [SVector{N,T}(x) for x in Codomain]
        else
            Codomain = [SVector{N,T}(Codomain)]
        end
    elseif Codomain isa Matrix{<:Number} # turn a matrix rowwise into a vector of SVector.
        Codomain = [SVector{N,T}(x) for x in eachrow(Codomain)]
    end

    _IsNormalized = CheckIfNormal(Domain, L) && CheckIfNormal(Codomain, L)

    return Crystal{N, T}(L, Domain, Codomain, _IsNormalized)
end

function Base.getproperty(C::Crystal, sym::Symbol)
    if sym == :size_domain
        length(C.Domain)
    elseif sym == :size_codomain
        length(C.Codomain)
    elseif sym ∈ [:A, :dim, :n, :iA, :dA] # some properties from lattice
        getproperty(getfield(C, :L), sym)
    else
        # fallback to getfield
        getfield(C, sym)
    end
end

"""
    normalize(C::Crystal{N,T}) where {N,T}

Normalizes the crystal, i.e., shifts the structure elements Domain and Codomain into the standard primitive cell.

# Example
```jldoctest
julia> using ALFA

julia> C = ALFA.Crystal{1,Float64}([2], [[-1],[-.5],[-1.5]], [[0]])
Lattice Basis: ALFA.Lattice{1,Float64}([2.0])
Domain: 3-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [-1.0]
 [-0.5]
 [-1.5]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]

julia> ALFA.normalize(C)
Lattice Basis: ALFA.Lattice{1,Float64}([2.0])
Domain: 3-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.5]
 [1.0]
 [1.5]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{1},Float64,1,1},1}:
 [0.0]

```
"""
function normalize(C::Crystal{N,T}) where {N,T}
    dn, ds, dp = ShiftIntoStandardCell(C.Domain, C.L)
    cn, ds, dp = ShiftIntoStandardCell(C.Codomain, C.L)
    return Crystal{N,T}(C.L, dn, cn)
end

"""
    wrtLattice(C::Crystal, L::Lattice)
    wrtLattice(C::Crystal{N,T}, A) where {N,T}

Rewrites the crystal with respect to L::Lattice or lattice basis A. Thus, the lattice
must be a sublattice of C.L.

# Example
```jldoctest
julia> using ALFA

julia> using LinearAlgebra

julia> C = ALFA.Crystal{2,Float64}()
Lattice Basis: ALFA.Lattice{2,Float64}([1.0 0.0; 0.0 1.0])
Domain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
Codomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]

julia> L = ALFA.Lattice{2,Float64}(2*I)
ALFA.Lattice{2,Float64}([2.0 0.0; 0.0 2.0])

julia> ALFA.wrtLattice(C,L)
Lattice Basis: ALFA.Lattice{2,Float64}([2.0 0.0; 0.0 2.0])
Domain: 4-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
 [1.0, 0.0]
 [0.0, 1.0]
 [1.0, 1.0]
Codomain: 4-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:
 [0.0, 0.0]
 [1.0, 0.0]
 [0.0, 1.0]
 [1.0, 1.0]

```
"""
function wrtLattice(C::Crystal{N,T}, A) where {N,T}
    t = ElementsInQuotientSpace(C.A, A, return_fractional = false)
    newDomain = [x + y for x in t for y in C.Domain]
    newCodomain = [x + y for x in t for y in C.Codomain]
    return Crystal{N,T}(A, newDomain, newCodomain)
end

function wrtLattice(C::Crystal, L::Lattice)
    return wrtLattice(C, L.A)
end

function Base.:(==)(A::Crystal, B::Crystal)
     return A.L == B.L && A.Domain == B.Domain && A.Codomain == B.Codomain
end

function Base.:(≈)(A::Crystal, B::Crystal)
     return A.L ≈ B.L && A.Domain ≈ B.Domain && A.Codomain ≈ B.Codomain
end
