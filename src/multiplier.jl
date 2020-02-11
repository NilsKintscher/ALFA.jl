struct Multiplier #<: AbstractArray{Int,1}
    pos::Vector{Int}
    mat::Matrix
    function Multiplier(pos::Vector{Int}, mat::Matrix{T}) where {T<:Number}
        new(pos, mat)
    end
end

Base.isless(a::Multiplier, b::Multiplier) = Base.isless(a.pos, b.pos)
Base.isequal(a::Multiplier, b::Multiplier) = Base.isequal(a.pos, b.pos)

function Multiplier(pos = nothing, mat = nothing)
    if pos == nothing
        pos = [0]
    elseif typeof(pos) == Matrix{Int} && size(pos, 1) == 1
        pos = reshape(pos, size(pos, 2)) # turning matrix with 1row into vector.
    end
    if mat == nothing
        mat = Matrix{Complex}(I, 0, 0)
    elseif typeof(mat) <: Vector
        mat = reshape(mat, length(mat), 1)
    end
    return Multiplier(convert(Array{Int,1}, pos), convert(Matrix, mat))
end

function Base.getproperty(m::Multiplier, sym::Symbol)
    if sym == :n || sym == :dim
        length(m.pos)
    else
        # fallback to getfield
        getfield(m, sym)
    end
end

function Base.size(m::Multiplier,y...)
    return size(m.mat,y...)
end
