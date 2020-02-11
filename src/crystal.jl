struct Crystal
    L::Lattice
    Domain::Matrix{Int} # m × dim
    Codomain::Matrix{Int} # m × dim
    function Crystal(L::Lattice, Domain::Matrix{Int}, Codomain::Matrix{Int})
        @assert L.dim == size(Domain, 2) "size(Domain,2)=$(size(Domain,2)) must be equal to Lattice dimensionality L.dim=$(L.dim)"
        @assert L.dim == size(Codomain, 2) "size(Codomain,2)=$(size(Codomain,2)) must be equal to Lattice dimensionality L.dim=$(L.dim)"
        if pointer(Domain) == pointer(Codomain)
             ## Or maybe dont do that and exploit it. Don't know yet.
            Codomain = deepcopy(Codomain)
        end
        new(L, Domain, Codomain)
    end
end

function Crystal(L = nothing, Domain = nothing, Codomain = nothing)
    if L == nothing
        L = Lattice()
    elseif typeof(L) != Lattice
        L = Lattice(L)
    end

    if Domain == nothing
        Domain = zeros(Int, 1, L.dim)
    elseif typeof(Domain) == Vector{Int} # turn Vector into Matrix with 1 Column
        Domain = reshape(Domain, length(Domain), 1)
    end

    if Codomain == nothing
        Codomain = Domain
    elseif typeof(Codomain) == Vector{Int} # turn Vector into Matrix with 1 Column
        Codomain = reshape(Codomain, length(Codomain), 1)
    end

    return Crystal(
        L,
        convert(Matrix{Int}, Domain),
        convert(Matrix{Int}, Codomain),
    )
end
