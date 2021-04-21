struct CrystalTorus{N,T}
    C::Crystal{N, T} # only using domain=codomain.
    Z::Lattice{N, T} # defines the torus size.
    coords::Vector{SVector{N,T}} # M::SortedSet{Multiplier}# Array{Multiplier,1}
    function CrystalTorus{N,T}(
        C::Crystal{N,T},
        Z::Lattice{N, T},
    ) where {N,T<:Union{Float64,Rational}}
        coords = ElementsInQuotientSpace(
            C.L.A,
            Z.A
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
