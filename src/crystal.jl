struct Crystal{N}
    L::Lattice{N}
    Domain::Vector{SVector{N,Float64}}   #Matrix # m × dim
    Codomain::Vector{SVector{N,Float64}} # m × dim
    _IsNormalized::Bool
    function Crystal{N}(
        L::Lattice{N},
        Domain::Vector{SVector{N,Float64}},
        Codomain::Vector{SVector{N,Float64}},
        _IsNormalized::Bool
    ) where {N}#
        new{N}(L, Domain, Codomain, _IsNormalized)
    end
end

function Crystal(L = nothing, Domain = nothing, Codomain = nothing)
    if L == nothing
        L = Lattice()
    elseif !isa(L, Lattice) #typeof(L) != Lattice
        L = Lattice(L)
    end

    N = typeof(L).parameters[1] # dimension

    if Domain == nothing # put a point at the origin.
        Domain = [zeros(SVector{N,Float64})]
    elseif typeof(Domain) <: Vector{Vector{T}} where {T<:Real} ||
           Domain isa Vector{Matrix{T}} where {T<:Real} #  turn Vector of Vector into Vector of SVector
        Domain = [SVector{N,Float64}(x) for x in Domain]
    elseif Domain isa Vector{<:Number} # turn vector into vector of SVector.
        if N == 1
            Domain = [SVector{N,Float64}(x) for x in Domain]
        else
            Domain = [SVector{N,Float64}(Domain)]
        end
    elseif Domain isa Matrix{<:Number} # turn a matrix rowwise into a vector of SVector.
        Domain = [SVector{N,Float64}(x) for x in eachrow(Domain)]
    end

    if Codomain == nothing # put a point at the origin.
        Codomain = Domain
    elseif typeof(Codomain) <: Vector{Vector{T}} where {T<:Real} ||
           Codomain isa Vector{Matrix{T}} where {T<:Real}# turn Vector into Matrix with 1 Column
        Codomain = [SVector{N,Float64}(x) for x in Codomain]
    elseif Codomain isa Vector{<:Number} # turn vector into vector of SVector.
        if N == 1
            Codomain = [SVector{N,Float64}(x) for x in Codomain]
        else
            Codomain = [SVector{N,Float64}(Codomain)]
        end
    elseif Codomain isa Matrix{<:Number} # turn a matrix rowwise into a vector of SVector.
        Codomain = [SVector{N,Float64}(x) for x in eachrow(Codomain)]
    end

    _IsNormalized = CheckIfNormal(Domain, L) && CheckIfNormal(Codomain, L)

    return Crystal{N}(L, Domain, Codomain, _IsNormalized)
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


function Base.show(io::IO, mime::MIME"text/plain", C::Crystal)
    print(io, "Lattice Basis: ")
    show(io, mime, C.L)
    print(io, "\nDomain: ")
    show(io, mime, C.Domain)
    print(io, "\nCodomain: ")
    show(io, mime, C.Codomain)
end


function normalize(C::Crystal)
    dn, ds, dp = ShiftIntoUnitCell(C.Domain, C.L)
    cn, ds, dp = ShiftIntoUnitCell(C.Codomain, C.L)
    return Crystal(C.L, dn, cn)
end

function wrtLattice(C::Crystal, A) # A::Matrix
    t = ElementsInQuotientSpace(C.A, A, return_fractional = false)
    newDomain = [x + y for x in t for y in C.Domain]
    newCodomain = [x + y for x in t for y in C.Codomain]
    return Crystal(A, newDomain, newCodomain)
end

function wrtLattice(C::Crystal, L::Lattice)
    return wrtLattice(C, L.A)
end

#
