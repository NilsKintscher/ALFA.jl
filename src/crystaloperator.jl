struct CrystalOperator
    C::Crystal ## dimension L.dim dictates size of Multiplier.pos and dim of structure elements.
    M::SortedSet{Multiplier}# Array{Multiplier,1}
    function CrystalOperator(C::Crystal, M::SortedSet{Multiplier})
        new(C,M)
    end
end

function CrystalOperator(C = nothing, M = nothing)
    if C == nothing
        C = Crystal()
    end
    if M == nothing
        M = SortedSet(Vector{Multiplier}())
    end
    return CrystalOperator(C, M)
end

function Base.push!(S::CrystalOperator, m::Multiplier)
    # check if already existent -> warning;  and check size: pos & matrix
    push!(S.M, m)
    return S
end
