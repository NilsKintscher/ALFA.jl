mutable struct CrystalVector{N,T,outerdim, innerdim}
    CT::CrystalTorus{N,T}
    v::MVector{outerdim, MVector{innerdim, Float64}}
    function CrystalVector{N,T, outerdim, innerdim}(
        CT::CrystalTorus{N,T},
        v::MVector{outerdim, MVector{innerdim, Float64}},
    ) where {N,T,innerdim, outerdim} # ]{N,T<:Union{Float64,Rational}, innerdim <: Int, outerdim <: Int}
        new{N,T, outerdim, innerdim}(CT, v)
    end
end

function CrystalVector(
    CT::CrystalTorus{N,T},
    v
    ) where {N,T}

    outerdim = length(CT.coords)
    innerdim = length(CT.C.Domain) # Codomain is not used at all.

    if typeof(v) <: Union{Vector, MVector}
        mv = MVector{outerdim}([
        MVector{innerdim}(vi) for vi in v
        ])
    else
        mv = MVector{outerdim}([
        MVector{innerdim}(v(innerdim)) for i in 1:outerdim
        ])
    end
    #return mv
    CrystalVector{N,T,outerdim, innerdim}(CT, mv)
end



function wrtLattice(CV::CrystalVector{N,T, outerdim, innerdim}, A) where {N,T, outerdim, innerdim}

    CT = CV.CT

    t = ElementsInQuotientSpace(CT.C.A, A, return_fractional = false)

    # construct new Torus
    newDomain = [x + y for x in t for y in CT.C.Domain]
    newCodomain = [x + y for x in t for y in CT.C.Codomain]

    C = Crystal{N,T}(A, newDomain, newCodomain)

    CTnew = CrystalTorus{N,T}(C, CT.Z)

    # todo: construct vector

    #
    return 0
end
