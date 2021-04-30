mutable struct CrystalVector{N,T,outerdim, innerdim}
    CT::CrystalTorus{N,T}
    v::SizedMatrix{outerdim, innerdim}
    function CrystalVector{N,T, outerdim, innerdim}(
        CT::CrystalTorus{N,T},
        v::SizedMatrix{outerdim, innerdim},
    ) where {N,T,outerdim, innerdim}
        new{N,T, outerdim, innerdim}(CT, v)
    end
end

function Base.getproperty(CV::CrystalVector{N,T,outerdim, innerdim}, sym::Symbol) where {N,T,outerdim, innerdim}
    if sym == :parameters
        (N,T,outerdim, innerdim)
    else
        # fallback to getfield
        getfield(CV, sym)
    end
end


function CrystalVector(
    CT::CrystalTorus{N,T},
    v
    ) where {N,T}

    outerdim = length(CT.coords)
    innerdim = length(CT.C.Domain) # Codomain is not used at all.

    if typeof(v) <: Union{Matrix, SizedMatrix}
        mv = SizedMatrix{outerdim, innerdim}(v)
        #SizedVector{outerdim}([
        #SizedVector{innerdim}(vi) for vi in v
        #])
    else
        mv = SizedMatrix{outerdim, innerdim}(v(outerdim, innerdim))
        #SizedVector{outerdim}([
        #SizedVector{innerdim}(v(innerdim)) for i in 1:outerdim
        #])
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
    # cnj are given in fractional coordinates in L(A)/L(Z)
    # the index j fulfils cnj - x_1 in L(Z)

    # construct new vector
    innerdim_new = innerdim*length(t)
    outerdim_new = length(CTnew.coords)
    v_new = SizedMatrix{outerdim_new, innerdim_new, typeof(CV.v).parameters[2]}(undef)
    #v_new[1,1] = 0

    #@show typeof(v_new)
     # typeof(CV.v).parameters[2]
    #SizedVector{outerdim_new,SizedVector{innerdim_new}}([SizedVector{innerdim_new}(undef) for i in 1:outerdim_new])

    #@show typeof(v_new)

    bs = CT.C.size_domain # blocksize

    for (it_c, c) in enumerate(CT.coords)
        #@show "NEW COORD"
        for (it_ti,ti) in enumerate(t)
            y = A\(CT.C.A*c-ti) # c and ti are fractional coordinates. transforming them to cart. coord and check if they are element of L(A)
            yr = round.(y)
            if isapprox(y, yr, rtol = ALFA_rtol, atol = ALFA_atol)
                # find the corresponding coordinate.
                for (it_cnew, cnew) in enumerate(CTnew.coords)
                    # cnew is given in fractional coordinates. transform to cartesian coordinates
                    # given the cartesian coordinate of y, we check if the difference is part of L(Z)
                    x = CTnew.Z.A\(CTnew.C.A*cnew - A*yr)
                    xr = round.(x)
                    if isapprox(x, xr, rtol = ALFA_rtol, atol = ALFA_atol)
                        #@show "FOUND"
                        #@show CV.v[it_c,:]
                        #@show v_new[it_cnew, (it_ti-1)*bs+1:it_ti*bs]
                        #v_new[1,1] = 0
                        #for i in 1:bs
                    #        v_new[it_cnew, (it_ti-1)*bs+i ] = 1.0 # CV.v[it_c,i] #1:it_ti*bs]
                #        end
                        #@show it_cnew, (it_ti-1)*bs+1:it_ti*bs
                        v_new[it_cnew, (it_ti-1)*bs+1:it_ti*bs] = CV.v[it_c,:]
                        #@show CV.v[it_c,:]
                        #@show v_new[it_cnew, (it_ti-1)*bs+1:it_ti*bs]
                        break
                    end
                end
                break
            end
        end
    end
    #@show typeof(v_new)
    #@show v_new[1]
    #return v_new
    #@show v_new
    CrystalVector{N,T,outerdim_new, innerdim_new}(CTnew, v_new)
    #println("a")
    #return 0
end


function Base.:(≈)(AV::CrystalVector{N,T,outerdim, innerdim}, BV::CrystalVector{N,T,outerdim, innerdim}) where {N,T,outerdim, innerdim}
    return AV.CT ≈ BV.CT && AV.v ≈ BV.v
end

function IsApproxEquivalent(A::CrystalVector, B::CrystalVector)
    (Anew, Bnew) = wrtSameLatticeAndNormalize(A, B)
    return Anew ≈ Bnew
end

function wrtSameLatticeAndNormalize(A::CrystalOperator, B::CrystalVector)
    B, A = wrtSameLatticeAndNormalize(B, A)
    return A, B
end
function wrtSameLatticeAndNormalize(AV::CrystalVector, S::CrystalOperator)
    if AV.CT.C.L.A == S.C.L.A
        Anew = AV
        Bnew = S
    else
        X = lcm(AV.CT.C.L, S.C.L)

        Anew = wrtLattice(AV, X)
        Bnew = wrtLattice(S, X)
    end
    Anew = normalize(Anew)
    Bnew = normalize(Bnew)
    return Anew, Bnew
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
        # dn_to in in unit cell, ds_to = shift.
        dn = Domain
        dp = dp_from[invperm(dp_to)]
        ds = ds_from-ds_to[invperm(dp_to)]

        for j in 1:length(Domain)
            @assert dn[j] + S.C.L.A*ds[j] ≈ S.C.Domain[dp[j]] "Domain structure element not compatible?, j= $j"
        end


        Cnew = Crystal{N,T}(S.C.L.A, dn)
        CTnew = CrystalTorus{N,T}(Cnew, CV.CT.Z, CV.CT.coords)

        v_new = deepcopy(CV.v)

        rowperm(x) = [mod(y-1+(x-1), outerdim)+1 for y in 1:outerdim]
        #@show rowperm(1)
        for (it,x) in enumerate(ds_to)
            # find CT
            for (it_c, c) in enumerate(CV.CT.coords)
                y = CV.CT.Z.A\(CV.CT.C.L.A*(x-c))
                yr = round.(y)
                if isapprox(y, yr, rtol = ALFA_rtol, atol = ALFA_atol)
                    #@show it
                    #@show it_c
                    #@show rowperm(it_c)
                    #@show v_new[:, it]
                    #@show v_new[:, it]
                    #@show "permute..:"
                    v_new[:, it] = v_new[rowperm(it_c), it]
                    #@show v_new[:, it]
                end
            end
        end
        v_new[:,:] = v_new[:,dp]
        CrystalVector{N,T,outerdim,innerdim}(CTnew, v_new)
    end



# Needed:
# Given a torus L(A)/L(Z) with coords:
# Shift coords into Z[0,1)^n

function ShiftCoordsIntoStandardCell(CV::CrystalVector{N,T,outerdim, innerdim}) where {N,T,outerdim, innerdim}
    s = [CV.CT.C.A*x for x in CV.CT.coords]

    t, y, p = ALFA.ShiftIntoStandardCell(s, CV.CT.Z)

    # t[j] + Z*y[j] = s[p[j]]
    # t[j] is found in the Standard Cell, cartesian coordinates.

    coords_new = [round.(CV.CT.C.A\x) for x in t]
    #coords_new = round.(CV.CT.C.A.\t) #fractional coordinates.

    v_new = CV.v[p,:] # SizedVector{outerdim}(CV.v[p])

    CT_new = CrystalTorus(CV.CT.C, CV.CT.Z, coords_new)

    #@show typeof(v_new)
    #return v_new
    CrystalVector(CT_new, v_new)
end
#SizedVector{9, SizedVector{4, TData} where TData<:AbstractVector{Float64}, Vector{SizedVector{4, TData} where TData<:AbstractVector{Float64}}}


function ChangeTorusCoords(CV::CrystalVector{N,T}, coords; fractional=true) where {N,T}
    # coords given in fractional coordinates?
    if(fractional)
        s_new_cart = [CV.CT.C.A*x for x in coords]
        s_new_fract = coords
    else
        s_new_cart = coords
        s_new_fract = [round.(CV.CT.C.A\x) for x in coords]
    end

    s_old_cart = [CV.CT.C.A*x for x in CV.CT.coords]

    (_, _, p_from) = ShiftIntoStandardCell(s_old_cart, CV.CT.Z)
    (_, _, p_to) = ShiftIntoStandardCell(s_new_cart, CV.CT.Z)


    v_new = CV.v[p_from[invperm(p_to)], :]
    CT_new = CrystalTorus(CV.CT.C, CV.CT.Z, s_new_fract)
    CrystalVector(CT_new, v_new)
end


# mathematical operations
function LinearAlgebra.norm(A::CrystalVector{N,T,outerdim, innerdim}) where {N,T,outerdim, innerdim}
    norm(A.v)
    #sum = 0
    #for vi in A.v
    #    sum += norm(vi)
    #end
    #sum/outerdim
end

# standalone
function Base.:-(A::CrystalVector)
    Ac = deepcopy(A)
    Ac.v = -Ac.v
    return Ac
end

# scalar:

function Base.:/(A::CrystalVector, b::T) where {T<:Number}
    return A * (1 / b)
end

function Base.:*(A::CrystalVector, b::T) where {T<:Number}
    Ac = deepcopy(A)
    Ac.v = Ac.v * b
    return Ac
end

function Base.:*(b::T, A::CrystalVector) where {T<:Number}
    return A*b
end

function Base.:+(b::T, A::CrystalVector) where {T<:Number}
    return A+b
end


function Base.:-(b::T, A::CrystalVector) where {T<:Number}
    return b+(-A)
end



function Base.:-(A::CrystalVector, b::T) where {T<:Number}
    return A+(-b)
end


function Base.:+(A::CrystalVector, b::T) where {T<:Number}
    Ac = deepcopy(A)
    Ac.v .+= b # = SizedVector([x.+b for x in Ac.v])
    return Ac
end

# vector with vector

function Base.:+(AV::CrystalVector, BV::CrystalVector)
    AV1, BV1 = wrtSameLatticeAndNormalize(AV,BV)
    @assert AV1.CT ≈ BV1.CT "Crystal Tori not equivalent."

    vnew = AV1.v + BV1.v
    CrystalVector(AV.CT, vnew)
end

function Base.:-(AV::CrystalVector, BV::CrystalVector)
    AV+(-BV)
end

function LinearAlgebra.dot(AV::CrystalVector, BV::CrystalVector)
    AV1, BV1 = wrtSameLatticeAndNormalize(AV,BV)
    @assert AV1.CT ≈ BV1.CT "Crystal Tori not equivalent."

    dot(AV1.v, BV1.v)
end

# CrystalOperator with vector
function Base.:*(A::CrystalOperator{N,T}, CV::CrystalVector{N,T,outerdim, innerdim}) where {N,T, outerdim, innerdim}
    # compute A*f = sum_{y in A.M} y.mat * f(x + y.pos)
    # todo
    A, CV = wrtSameLatticeAndNormalize(A,CV)
    @assert A.C.L ≈ CV.CT.C.L "Lattices not equivalent."
    @assert A.C.Domain ≈ CV.CT.C.Domain "Domains not equivalent."
    #
    Cnew = Crystal{N,T}(CV.CT.C.L, A.C.Codomain)
    CTnew = CrystalTorus(Cnew,CV.CT.Z, CV.CT.coords)

    outerdim_new = CV.parameters[3] # outerdim
    innerdim_new = length(A.C.Codomain)


    #v_new = SizedVector{outerdim_new,SizedVector{innerdim_new}}([SizedVector{innerdim_new}(zeros(innerdim_new)) for i in 1:outerdim_new])
    v_new = SizedMatrix{outerdim_new, innerdim_new}(zeros(outerdim_new, innerdim_new))

    for (v_pos_frac_it, v_pos_frac) in  enumerate(CV.CT.coords)
        #vi = zeros(innerdim_new)
        for am in A.M
            # m_pos = am.pos
            # m_mat = am.mat

            # find x = v_pos_frac + y.pos in coords.
            # x is given in fractional coordinates.
            pos_frac = v_pos_frac + am.pos
            #@show pos_frac
            found = false
            for (v_from_pos_frac_it, v_from_pos_frac) in enumerate(CV.CT.coords)
                x = CV.CT.Z.A\(A.C.L.A*(pos_frac - v_from_pos_frac))
                xr = round.(x)
                if isapprox(x, xr, rtol = ALFA_rtol, atol = ALFA_atol)
                    #@show CV.v[v_from_pos_frac_it]
                    #@show am.mat
                    #@show typeof(am.mat*CV.v[v_from_pos_frac_it])
                    #vi += am.mat*CV.v[v_from_pos_frac_it]
                    v_new[v_pos_frac_it,:] += am.mat*CV.v[v_from_pos_frac_it,:]
                    found = true
                    break
                end
            end
            @assert found "No corresponding position found. Something must be wrong."
        end
    end
    #return CTnew, v_new
    CrystalVector(CTnew, v_new)
end


#
