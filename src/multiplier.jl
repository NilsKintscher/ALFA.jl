mutable struct Multiplier{N} #<: AbstractArray{Int,1}
    pos::MVector{N}
    mat::Matrix
    function Multiplier{N}(pos::MVector{N,Int}, mat::Matrix{T}) where {N,T}
        new{N}(pos, mat)
    end
end

Base.isless(a::Multiplier, b::Multiplier) = Base.isless(a.pos, b.pos)
Base.isequal(a::Multiplier, b::Multiplier) = Base.isequal(a.pos, b.pos)

"""
    Multiplier{N}(pos = nothing, mat = nothing) where {N}
    Multiplier(pos = nothing, mat = nothing)

Constructs a multiplication matrix as part of a CrystalOperator.
The position is given in fractional coordinates and thus (converted to) integral.

# Example
```jldoctest
julia> using ALFA

julia> ALFA.Multiplier([0 0], [1 2 3; 4 5 6])
Position: 2-element StaticArraysCore.MVector{2, Int64} with indices SOneTo(2):
 0
 0
Multiplier: 2×3 Matrix{Int64}:
 1  2  3
 4  5  6

```
"""
function Multiplier{N}(pos = nothing, mat = nothing) where {N}
    return Multiplier(pos, mat)
end

function Multiplier(pos = nothing, mat = nothing)
    if pos === nothing
        pos = MVector{1,Int}(0)
    end
    N = length(pos)
    if mat === nothing
        mat = Matrix{Complex}(I, 0, 0)
    elseif typeof(mat) <: Vector
        mat = reshape(mat, length(mat), 1)
    end

    return Multiplier{N}(convert(MVector{N,Int}, pos), convert(Matrix, mat))
end

function Base.getproperty(m::Multiplier, sym::Symbol)
    if sym == :n || sym == :dim
        length(m.pos)
    elseif sym == :size_domain
        size(m.mat, 2)
    elseif sym == :size_codomain
        size(m.mat, 1)
    else
        # fallback to getfield
        getfield(m, sym)
    end
end

function Base.size(m::Multiplier, y...)
    return size(m.mat, y...)
end

function Base.:(==)(A::Multiplier, B::Multiplier)
     return A.pos == B.pos && A.mat == B.mat
end

function Base.:(≈)(A::Multiplier, B::Multiplier)
     return A.pos ≈ B.pos && A.mat ≈ B.mat
end
