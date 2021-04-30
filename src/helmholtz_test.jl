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



function all_combs(CV::CrystalVector{N,T, outerdim, innerdim}, n::Int) where {N,T, outerdim, innerdim}

CT4 = wrtLattice(CT, CT.C.L.A*n)

xZ4 = ALFA.wrtCrystalReverse(xZ, CT4)
xZ4_all = [ALFA.ChangeStructureElement(xZ4, [y .+ x for y in CT4.C.Domain])  for x in CT4.C.Domain ]

end
