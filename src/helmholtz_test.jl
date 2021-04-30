#
function wrtCrystalReverse(CV::CrystalVector{N,T, outerdim, innerdim}, CT::CrystalTorus{N,T}) where {N,T, outerdim, innerdim}

    CV = ALFA.wrtLattice(CV, CV.CT.Z) # unwrap vector
    CV = ALFA.normalize(CV)

    CT = ALFA.normalize(CT)

    innerdim_new = length(CT.C.Domain)
    outerdim_new = length(CT.coords)

    v_new = SizedMatrix{outerdim_new, innerdim_new, typeof(CV.v).parameters[2]}(undef)

    for (c_it, c) in enumerate(CT.coords)
        c_cart = CT.C.L.A*c
        for (d_it, d) in enumerate(CT.C.Domain)
            x = c_cart + d # corresponds to some position in old domain modulo Z
            for (c_from_it, c_from) in enumerate(CV.CT.C.Domain)
                y = CV.CT.Z.A\(x - c_from)
                yr = round.(y)
                if isapprox(y, yr, rtol = ALFA_rtol, atol = ALFA_atol)
                    v_new[c_it, d_it] = CV.v[1, c_from_it]
                    break
                end
            end
        end
    end
    CrystalVector(CT, v_new)
end



function all_combs(CV::CrystalVector{N,T, outerdim, innerdim}, CT::CrystalTorus{N,T}) where {N,T, outerdim, innerdim}

xZ = ALFA.wrtCrystalReverse(CV, CT)

xZ_all = [ALFA.ChangeStructureElement(xZ, [y .+ x for y in CT.C.Domain])  for x in CT.C.Domain ]


# reinterpret them wrt to the original structure element.
return [ALFA.CrystalVector(xZ.CT, x.v) for x in xZ_all]
end

function svd_coarsening(CVfvec,num_coarse_points)
    #assuming that the underlying structures of CVfvec are identical.
    val_distr_fine = vcat([x.v for x in CVfvec]...)
    U,S,V = svd(val_distr_fine)

    val_distr_coarse = val_distr_fine * V[:,1:num_coarse_points]

    # construct coarse vectors
    Cf = CVfvec[1].CT.C
    A = Cf.L.A


    Domain = Cf.Domain[1:Int(floor(length(Cf.Domain)/num_coarse_points)):end]

    C = ALFA.Crystal{Cf.n,Float64}(A, Domain)

    CVcvec = []
    bs = length(CVfvec[1].CT.coords)
    for (cvf_it,cvf_i) in enumerate(CVfvec)
        CT = ALFA.CrystalTorus(C, CVfvec[1].CT.Z)
        idx_begin = (cvf_it-1)*bs+1
        idx_end = (cvf_it)*bs

        v = val_distr_coarse[idx_begin:idx_end, : ]
        push!(CVcvec, ALFA.CrystalVector(CT, v))
    end
    return CVcvec
end


@with_kw struct interp_strategy_params
           idx_codomain_vec = [[1,2], [3,4]]
           idx_domain_vec = [[1,2], [1,2]]
           coord_vec_list = [[[-1], [0]], [[0], [1]]]
       end


_idx_1d_4_to_2 = interp_strategy_params([[1,2], [3,4]],
    [[1,2], [1,2]],
    [[[-1], [0]], [[0], [1]]]
    )


function construct_P(CVfvec, CVcvec, interp_stategy::interp_strategy_params=_idx_1d_4_to_2)
    A = CVfvec[1].CT.C.L.A
    CT = CVfvec[1].CT
    Codomain = CVfvec[1].CT.C.Domain
    Domain = CVcvec[1].CT.C.Domain

    C = ALFA.Crystal{CVfvec[1].CT.C.n, Float64}(A, Domain, Codomain)

    P  = ALFA.CrystalOperator{CVfvec[1].CT.C.n,Float64}(C)

    # what idx to prolongate to (idx of fine vector.)
    #idx_codomain_vec = [[1,2], [3,4]]
    #idx_codomain = [1,2]
    # from where (idx and corr. pos of coarse vector)
    #idx_domain_vec = [[1,2], [1,2]]
    #idx_domain = [1,2]
    #coord_vec_list = [[[-1], [0]], [[0], [1]]]
    #coord_vec = [[-1], [0]] # fractional coordinates
    #@show "1"
    #idx_domain_vec, idx_codomain_vec, coord_vec_list = interp_strategy()
    # interp_strategy

    #idx_domain_vec = _idx_1d_4_to_2.idx_domain_vec
    #idx_codomain_vec = _idx_1d_4_to_2.idx_codomain_vec
    #coord_vec_list = _idx_1d_4_to_2.coord_vec_list
    @unpack idx_domain_vec, idx_codomain_vec, coord_vec_list = interp_stategy


    for it in 1:length(idx_codomain_vec)
        idx_codomain = idx_codomain_vec[it]
        idx_domain = idx_domain_vec[it]
        coord_vec = coord_vec_list[it]
        weights = _get_weights(idx_domain, idx_codomain, coord_vec, CT, CVfvec, CVcvec)


    bs = length(idx_domain)
    for (pos_it, pos) in enumerate(coord_vec)
            m = zeros(C.size_codomain, C.size_domain)
            idx_begin = (pos_it-1)*bs+1
            idx_end = (pos_it)*bs
            m[idx_codomain,:] = weights[idx_begin:idx_end, :]'
            push!(P, ALFA.Multiplier(pos, m), true)
    end
end
    return P

end

function _get_weights( idx_domain, idx_codomain, coord_pos, CT, CVfvec, CVcvec)
    idx_vc = zeros(Int,length(CT.coords), length(coord_pos))
    # get for any coord the left neighbor
    #all vectors are normalized and describe the value distributions in the same way. Using the first vector as a representant.
    for (pos_frac_it, pos_frac) in  enumerate(CT.coords)
        for (shift_pos_it, shift_pos) in enumerate(coord_pos)
            x = pos_frac + shift_pos
            # find idx in coords
            for (ys_it, ys) in enumerate(CT.coords)
                y = CT.Z.A\(CT.C.L.A*(x-ys)) # to
                yr = round.(y)
                if isapprox(y, yr, rtol = ALFA_rtol, atol = ALFA_atol)
                    idx_vc[pos_frac_it, shift_pos_it] = ys_it
                    break
                end
            end
        end
    end
    vc = vcat([hcat( [x.v[idx_vc[:,i],:] for i in 1:length(coord_pos)]...) for x in CVcvec]...) # enumerate(vc_all)]
    vf = vcat([x.v[:,idx_codomain] for x in CVfvec]...)

    F = qr(vc)
    weights = F\vf
end
