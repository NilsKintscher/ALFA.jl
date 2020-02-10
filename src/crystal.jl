struct Crystal
    L::Lattice
    Domain::Vector{Vector{Int}}
    Codomain::Vector{Vector{Int}}
    function Crystal(L::Lattice, Domain::Vector{Vector{Int}}, Codomain::Vector{Vector{Int}})
        if pointer(Domain) == pointer(Codomain)
            Codomain = deepcopy(Codomain) ## Or maybe dont do that and exploit it. Don't know yet.
        end
        new(L,Domain,Codomain)
    end
end

function Crystal(L = nothing, Domain = nothing, Codomain = nothing)
    if L == nothing
        L = Lattice()
    elseif typeof(L) != Lattice
        L = Lattice(L) # try to call the constructor of Lattice
    end
    if Domain == nothing
        Domain = [[0, 0]]
    end
    if Codomain == nothing
        Codomain = Domain
    end
    return Crystal(L, Domain, Codomain)
end
