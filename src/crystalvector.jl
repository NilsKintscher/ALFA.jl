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
        MVector{innerdim, Float64}(vi) for vi in v
        ])
    else
        mv = MVector{outerdim}([
        MVector{innerdim,Float64}(v(innerdim)) for i in 1:outerdim
        ])
    end
    #return mv
    CrystalVector{N,T,outerdim, innerdim}(CT, mv)
end



function wrtLattice(CV::CrystalVector{N,T, outerdim, innerdim}, A) where {N,T, outerdim, innerdim}

    if typeof(A) <: ALFA.Lattice
        A = A.A
    end
    CT = CV.CT

    t = ElementsInQuotientSpace(CT.C.A, A, return_fractional = false)

    # construct new Torus
    newDomain = [x + y for x in t for y in CT.C.Domain]
    newCodomain = [x + y for x in t for y in CT.C.Codomain]

    C = Crystal{N,T}(A, newDomain, newCodomain)

    CTnew = CrystalTorus(C, CT.Z)

    # construct new vector

    #  t: coordinates in L(CT.C.A)/L(A) that will be added to the new structure element

    #  Need to find the CT.coords (which are part of L(CT.C.A)) that correspond to to one new structure element, i.e.
    # For each new Lattice Point in L(A)/L(CT.Z), there exist exactly Length(t) number of  CT.coords.
    # In other words:
    # CT.coords is divided in length(t) groups.
    # Each group corresponds to a value distribution of the structure element corresponding to exactly one lattice point of L(A)/L(CT.Z)


    # We first need to find the correspondance of c = CT.coords to t
    # for each c[i] exists exactly one t[j], such that the (c[i] - t[j]) in L(A)

    # example: CT.C.A = [1]
    # A = [2]
    # CT.Z = [4]
    # then
    # c = [0,1,2,3] coordinates in L(CT.C.A) / L(Z)
    # t = [0,1] coordinates in L(CT.C.A) / L(A)
    # and
    # For c[1]=0 we find c[1] - t[1] =x_1 in L(A)
    # For c[2]=1 we find c[2] - t[2] =x_2 in L(A)
    #
    # Thus, v[1] and v[2] are merged into one new vector vnew[j]=[c[1], c[2]] which corresponds to a single new coord cnj = CTnew.coords[j].
    # cnj are fractional coordinates in L(A)/L(Z)
    # the index j fulfils cnj - x_1 in L(Z)

    # construct new vector
    innerdim_new = innerdim*length(t)
    outerdim_new = length(CTnew.coords)
    v_new = MVector{outerdim_new,MVector{innerdim_new,Float64}}([MVector{innerdim_new,Float64}(undef) for i in 1:outerdim_new])

    bs = CT.C.size_domain # blocksize

    for (it_c, c) in enumerate(CT.coords)
#println("")
#@show CT.C.A
#@show A
#@show A\CT.C.A
#@show t

        #println("new c:")


        #@show (it_c, c) # fractional coordinate
        for (it_ti,ti) in enumerate(t)
            #@show (it_ti,ti)
            y = A\(CT.C.A*c-ti) # c and ti are fractional coordinates. transforming them to cart. coord and check if they are element of L(A)
            yr = round.(y)
            #@show y
            if isapprox(y, yr, rtol = ALFA_rtol, atol = ALFA_atol)
                #println("found ti")
                # find the corresponding coordinate.
                for (it_cnew, cnew) in enumerate(CTnew.coords)
                    #@show (it_cnew, cnew)
                    # cnew is given in fractional coordinates. transform to cartesian coordinates
                    # given the cartesian coordinate of y, we check if the difference is part of L(Z)
                    x = CTnew.Z.A\(CTnew.C.A*cnew - A*yr)
                    xr = round.(x)
                    #@show x
                    if isapprox(x, xr, rtol = ALFA_rtol, atol = ALFA_atol)
                        #println("found also cnew")
                        # found it.
                        #@show it_cnew, (it_ti-1)*bs+1:it_ti*bs
                        #@show it_c
                        v_new[it_cnew][(it_ti-1)*bs+1:it_ti*bs] = CV.v[it_c]
                        #@show "NEXT"
                        #println("")
                        break
                    end
                end
                break
            end
        end
    end

    CrystalVector{N,T,outerdim_new, innerdim_new}(CTnew, v_new)
end


function Base.:(≈)(AV::CrystalVector{N,T,outerdim, innerdim}, BV::CrystalVector{N,T,outerdim, innerdim}) where {N,T,outerdim, innerdim}
    if AV.CT ≈ BV.CT
        for i in 1:outerdim
            if !(AV.v[i] ≈ BV.v[i])
                return false
            end
        end
        return true
    else
        return false
    end
end

function IsApproxEquivalent(A::CrystalVector, B::CrystalVector)
    (Anew, Bnew) = wrtSameLatticeAndNormalize(A, B)
    return Anew ≈ Bnew
end


function wrtSameLatticeAndNormalize(AV::CrystalVector, BV::CrystalVector)

    if AV.CT.C.L.A == BV.CT.C.L.A
        Anew = AV
        Bnew = BV
    else
        X = lcm(AV.CT.C.L, BV.CT.C.L)
        
        Anew = wrtLattice(AV, X)
        Bnew = wrtLattice(BV, X)
    end
    Anew = normalize(Anew)
    Bnew = normalize(Bnew)
    return Anew, Bnew
end

function normalize(CV::CrystalVector{N,T,outerdim, innerdim}) where {N,T,outerdim, innerdim}
    S = CV.CT
    if !(S.C._IsNormalized)
    (dn, ds, dp) = ShiftIntoStandardCell(S.C.Domain, S.C.L)

    CV = ChangeStructureElement(CV, dn)
    end

    ShiftCoordsIntoStandardCell(CV)
end



    function ChangeStructureElement(CV::CrystalVector{N,T,outerdim, innerdim}, Domain) where {N,T, outerdim, innerdim}
        S = CV.CT
        (dn_from, ds_from, dp_from) = ShiftIntoStandardCell(S.C.Domain, S.C.L)

        newC = ALFA.Crystal{N,T}(S.C.L, Domain)
        (dn_to, ds_to, dp_to) = ShiftIntoStandardCell(newC.Domain, S.C.L)

        dn = Domain
        dp = dp_from[invperm(dp_to)]
        ds = ds_from-ds_to[invperm(dp_to)]

        for j in 1:length(Domain)
            @assert dn[j] + S.C.L.A*ds[j] ≈ S.C.Domain[dp[j]] "Domain structure element not compatible?, j= $j"
        end


        Cnew = Crystal{N,T}(S.C.L.A, dn)
        CTnew = CrystalTorus{N,T}(Cnew, CV.CT.Z, CV.CT.coords)
        #dn[j] = CT.C.Domain[dp[j]]
        #Thus, we need vnew[i][j] = v[i][dp[j]]
        v_new = copy(CV.v)
        for v in v_new
            v = v[dp]
        end

        CrystalVector{N,T,outerdim,innerdim}(CTnew, v_new)
    end



# Needed:
# Given a torus L(A)/L(Z) with coords:
# Shift coords into Z[0,1)^n

function ShiftCoordsIntoStandardCell(CV::CrystalVector)
    s = [CV.CT.C.A*x for x in CV.CT.coords]

    t, y, p = ALFA.ShiftIntoStandardCell(s, CV.CT.Z)

    # t[j] + Z*y[j] = s[p[j]]
    # t[j] is found in the Standard Cell, cartesian coordinates.

    coords_new = [round.(CV.CT.C.A\x) for x in t]
    #coords_new = round.(CV.CT.C.A.\t) #fractional coordinates.

    v_new = CV.v[p]

    CT_new = CrystalTorus(CV.CT.C, CV.CT.Z, coords_new)

    CrystalVector(CT_new, v_new)
end
