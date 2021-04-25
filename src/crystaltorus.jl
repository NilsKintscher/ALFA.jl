struct CrystalTorus{N,T}
    C::Crystal{N, T} # only using domain=codomain.
    Z::Lattice{N, T} # defines the torus size.
    coords::Vector{SVector{N,T}} # in fractional coordinates.
    function CrystalTorus{N,T}(
        C::Crystal{N,T},
        Z::Lattice{N, T},
        coords::Vector{SVector{N,T}}
    ) where {N,T<:Union{Float64,Rational}}
        new{N,T}(C, Z, coords)
    end
end

function CrystalTorus(
    C::Crystal{N,T},
    Z::Lattice{N, T},
    coords::Vector{SVector{N,T}}
) where {N,T<:Union{Float64,Rational}}
    CrystalTorus{N,T}(C, Z, coords)
end


function CrystalTorus(
    C::Crystal{N,T},
    Z::Lattice{N, T}
) where {N,T<:Union{Float64,Rational}}
    coords = ElementsInQuotientSpace(
        C.L.A,
        Z.A,
        return_fractional = true
    )
    CrystalTorus(C,Z,Vector{SVector{N,T}}(coords))
end


function CrystalTorus(
    C::Crystal{N,T},
    z::Int) where {N,T}
    Z = ALFA.Lattice{N,T}(z*C.L.A)
    CrystalTorus(C,Z)
end

function CrystalTorus(
    C::Crystal{N,T},
    Z::X) where {N,T,X<:Union{Matrix, MMatrix}}
    Z = ALFA.Lattice{N,T}(Z)
    CrystalTorus(C,Z)
end

# function CrystalTorus(
#     C::Crystal{N,T},
#     Z::ALFA.Lattice{N,T}) where {N,T}
#     CrystalTorus(C,Z)
# end

function wrtLattice(CT::CrystalTorus{N,T}, A) where {N,T}

    t = ElementsInQuotientSpace(CT.C.A, A, return_fractional = false)
    newDomain = [x + y for x in t for y in CT.C.Domain]
    newCodomain = [x + y for x in t for y in CT.C.Codomain]

    C = Crystal{N,T}(A, newDomain, newCodomain)

    return CrystalTorus{N,T}(C, CT.Z)
end


function Base.:(≈)(A::CrystalTorus, B::CrystalTorus)
    if A.C ≈ B.C && A.Z ≈ B.Z && A.coords ≈ B.coords
        return true
    else
        return false
    end
end


function normalize(CT::CrystalTorus{N,T}) where {N,T}
    CTnC = normalize(CT.C)
    return ALFA.CrystalTorus(CTnC, CT.Z)
end
