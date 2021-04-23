struct CrystalTorus{N,T}
    C::Crystal{N, T} # only using domain=codomain.
    Z::Lattice{N, T} # defines the torus size.
    coords::Vector{SVector{N,T}} # in fractional coordinates.
    function CrystalTorus{N,T}(
        C::Crystal{N,T},
        Z::Lattice{N, T},
    ) where {N,T<:Union{Float64,Rational}}
        coords = ElementsInQuotientSpace(
            C.L.A,
            Z.A,
            return_fractional = true
        )
        new{N,T}(C, Z, coords)
    end
end

function CrystalTorus(
    C::Crystal{N,T},
    z::Int) where {N,T}
    Z = ALFA.Lattice{N,T}(z*C.L.A)
    CrystalTorus{N,T}(C,Z)
end

function wrtLattice(CT::CrystalTorus{N,T}, A) where {N,T}

    t = ElementsInQuotientSpace(CT.C.A, A, return_fractional = false)
    newDomain = [x + y for x in t for y in CT.C.Domain]
    newCodomain = [x + y for x in t for y in CT.C.Codomain]

    C = Crystal{N,T}(A, newDomain, newCodomain)

    return CrystalTorus{N,T}(C, CT.Z)
end
