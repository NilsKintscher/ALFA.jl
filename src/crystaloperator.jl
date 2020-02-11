struct CrystalOperator
    C::Crystal ## dimension L.dim dictates size of Multiplier.pos and dim of structure elements.
    M::SortedSet{Multiplier}# Array{Multiplier,1}
    function CrystalOperator(C::Crystal, M::SortedSet{Multiplier})
        _sanitycheck(C,M)
        new(C, M)
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
