struct Crystal
    L::Lattice
    Domain::Matrix # m × dim
    Codomain::Matrix # m × dim
    function Crystal(L::Lattice, Domain::Matrix, Codomain::Matrix)
        @assert L.dim == size(Domain, 2) "size(Domain,2)=$(size(Domain,2)) must be equal to Lattice dimensionality L.dim=$(L.dim)"
        @assert L.dim == size(Codomain, 2) "size(Codomain,2)=$(size(Codomain,2)) must be equal to Lattice dimensionality L.dim=$(L.dim)"
        @assert typeof(Domain) <: Matrix{<:Real} "Domain must be of type <: Matrix{<:Real}"
        @assert typeof(Codomain) <: Matrix{<:Real} "Codomain must be of type <: Matrix{<:Real}"
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
    elseif typeof(Domain) <: Vector # turn Vector into Matrix with 1 Column
        Domain = reshape(Domain, length(Domain), 1)
    end

    if Codomain == nothing
        Codomain = Domain
    elseif typeof(Codomain) <: Vector # turn Vector into Matrix with 1 Column
        Codomain = reshape(Codomain, length(Codomain), 1)
    end

    return Crystal(L, convert(Matrix, Domain), convert(Matrix, Codomain))
end


function Base.getproperty(C::Crystal, sym::Symbol)
    if sym == :size_domain
        size(C.Domain, 1)
    elseif sym == :size_codomain
        size(C.Codomain, 1)
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

function normalize!(C::Crystal)
    ShiftIntoUnitCell!(C.Domain, C.L)
    ShiftIntoUnitCell!(C.Codomain, C.L)
end

function normalize(C::Crystal)
    dn, ds, dp = ShiftIntoUnitCell(C.Domain, C.L)
    cn, ds, dp = ShiftIntoUnitCell(C.Codomain, C.L)
    return Crystal(C.L, dn, cn)
end

function wrtLattice(C::Crystal, A::Matrix)
    t = ElementsInQuotientSpace(C.A, A, fractional=false)
    newDomain = vcat(transpose([x + y for x in eachslice(t, dims = 1) for y in eachslice(
        C.Domain,
        dims = 1,
    )])...)
    newCodomain = vcat(transpose([x + y for x in eachslice(t, dims = 1) for y in eachslice(
        C.Codomain,
        dims = 1,
    )])...)
    return Crystal(A, newDomain, newCodomain)
end

function wrtLattice(C::Crystal, L::Lattice)
    return wrtLattice(C, L.A)
end

# function Plots.plot(C::Crystal; xmin=-2, xmax=2, draw_basis=true)
#     Plots.plot()
#     plot!(C, xmin=xmin, xmax=max, draw_basis=draw_basis)
# end
#
# function Plots.plot!(C::Crystal; xmin=-2, xmax=2, ymin=-2, ymax=2, draw_basis=true)
#     plot!(C.L, xmin=xmin, xmax=max, draw_basis=draw_basis) ## TODO: handle kwargs in a better,simpler way
#     pos_fractional = Iterators.product(Iterators.repeated(xmin:xmax,C.dim)...) # fractional positions
#     plot!(C, pos_fractional)
# end
#
# function Plots.plot!(C::Crystal, pos_fractional; xmin=-2, xmax=2, ymin=-2, ymax=2, draw_basis=true)
#     pos_cartesian = [C.L.A*[x...] for x in pos_frac]
#     for x in pos_cartesian
#         scatter!([x+y for y in C.Domain])
#     end
#
# end
