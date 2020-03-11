mutable struct Multiplier{N} #<: AbstractArray{Int,1}
    pos::MVector{N}
    mat::Matrix
    function Multiplier{N}(pos::MVector{N,Int}, mat::Matrix{T}) where {N,T<:Number}
        new{N}(pos, mat)
    end
end

Base.isless(a::Multiplier, b::Multiplier) = Base.isless(a.pos, b.pos)
Base.isequal(a::Multiplier, b::Multiplier) = Base.isequal(a.pos, b.pos)

function Multiplier(pos = nothing, mat = nothing)
    if pos == nothing
        pos = MVector{1,Int}(0)
    # elseif pos isa Matrix || pos isa Vector
    #     pos = SVector{length(pos),Int}(pos)
    end
    N = length(pos)
    if mat == nothing
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
