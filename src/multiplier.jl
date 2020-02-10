struct Multiplier #<: AbstractArray{Int,1}
    pos::Vector{Int}
    mat::Matrix{Complex}
    function Multiplier(pos::Vector{Int}, mat::Matrix{Complex})
        new(pos, mat)
    end
end

Base.isless(a::Multiplier, b::Multiplier) = Base.isless(a.pos, b.pos)
Base.isequal(a::Multiplier, b::Multiplier) = Base.isequal(a.pos, b.pos)

function Multiplier(pos = nothing, mat = nothing)
    if pos == nothing
        pos = [0]
    end
    if mat == nothing
        mat = Matrix{Complex}(I, 0, 0)
    end
    return Multiplier(convert(Array{Int,1}, pos), convert(Matrix{Complex}, mat))
end
