#
module gallery
using ..alfa
using LinearAlgebra

function Laplace2D()
    A = [1 0; 0 1]
    Domain = [0 0]
    Codomain = [0 0]

    C = alfa.Crystal(A, Domain, Codomain)
    L = alfa.CrystalOperator(C)

    push!(L, alfa.Multiplier([0 0], [-4]))
    push!(L, alfa.Multiplier([0 -1], [1]))
    push!(L, alfa.Multiplier([0 1], [1]))
    push!(L, alfa.Multiplier([1 0], [1]))
    push!(L, alfa.Multiplier([-1 0], [1]))
    return L
end

function fw_restriction2D()
    A = 2 * [1 0; 0 1]
    Domain = [[0, 0], [0, 1], [1, 0], [1, 1]]
    Codomain = [0, 0]

    C = alfa.Crystal(A, Domain, Codomain)
    L = alfa.CrystalOperator(C)

    push!(L, alfa.Multiplier([0 0], [1 1 // 2 1 // 2 1 // 4]))
    push!(L, alfa.Multiplier([-1 0], [0 0 1 // 2 1 // 4]))
    push!(L, alfa.Multiplier([0 -1], [0 1 // 2 0 1 // 4]))
    push!(L, alfa.Multiplier([-1 -1], [0 0 0 1 // 4]))
    return L
end

function graphene_tight_binding(t = nothing)
    if t == nothing
        t = [0, -1, 0, 0]
    end
    A = 1 // 2 * [3 3; sqrt(3) -sqrt(3)]
    Domain = [j // 3 * A * [1, 1] for j in [1, 2]]

    C = alfa.Crystal(A, Domain)
    L = alfa.CrystalOperator(C)

    push!(L, alfa.Multiplier([0, 0], [t[1] t[2]; t[2] t[1]]))

    push!(L, alfa.Multiplier([-1, 0], [t[3] t[2]; 0 t[3]]))
    push!(L, alfa.Multiplier([0, -1], [t[3] t[2]; 0 t[3]]))

    push!(L, alfa.Multiplier([1, 0], [t[3] 0; t[2] t[3]]))
    push!(L, alfa.Multiplier([0, 1], [t[3] 0; t[2] t[3]]))

    push!(L, alfa.Multiplier([1, -1], [t[3] t[4]; t[4] t[3]]))
    push!(L, alfa.Multiplier([-1, 1], [t[3] t[4]; t[4] t[3]]))

    push!(L, alfa.Multiplier([1, 1], [0 0; t[4] 0]))
    push!(L, alfa.Multiplier([-1, -1], [0 t[4]; 0 0]))

    return L
end

function graphene_dirac_restriction(; m = 2, wl = 0.5, wlh = 0)
    ws = wl + wlh - 1

    A = 2^(m - 2) * 1 // 2 * [3 3; sqrt(3) -sqrt(3)]
    Ac = 2 * A

    s01 = [j // 3 * A * [1, 1] for j in [1, 2]]

    Domain = [A * x + y for x in [[0, 0], [1, 0], [0, 1], [1, 1]] for y in s01]
    Codomain = [j // 3 * Ac * [1, 1] for j in [1, 2]]

    C = alfa.Crystal(Ac, Domain, Codomain)
    L = alfa.CrystalOperator(C)

    m = zeros(2, 8)
    m[2, 1] = wl
    push!(L, alfa.Multiplier([1, 1], m))


    m = zeros(2, 8)
    m[2, [1, 3, 5]] = [ws, wlh, ws]
    push!(L, alfa.Multiplier([1, 0], m))

    m = zeros(2, 8)
    m[2, [1, 3, 5]] = [ws, ws, wlh]
    push!(L, alfa.Multiplier([0, 1], m))

    m = zeros(2, 8)
    m[2, 5] = wl
    m[1, 6] = wl
    push!(L, alfa.Multiplier([1, -1], m))

    m = zeros(2, 8)
    m[1, [2, 4, 6, 8]] = [1, ws, ws, wlh]
    m[2, [1, 3, 5, 7]] = [wlh, ws, ws, 1]
    push!(L, alfa.Multiplier([0, 0], m))

    m = zeros(2, 8)
    m[2, 3] = wl
    m[1, 4] = wl
    push!(L, alfa.Multiplier([-1, 1], m))

    m = zeros(2, 8)
    m[1, [4, 6, 8]] = [wlh, ws, ws]
    push!(L, alfa.Multiplier([0, -1], m))

    m = zeros(2, 8)
    m[1, [4, 6, 8]] = [ws, wlh, ws]
    push!(L, alfa.Multiplier([-1, 0], m))

    m = zeros(2, 8)
    m[1, 8] = wl
    push!(L, alfa.Multiplier([-1, -1], m))

    return L
end



function Base.rand(::Type{alfa.CrystalOperator}; single_domain = false, maxdigits = 3)
    MaxDim = 4 #4
    MaxNumElements = 2 #10
    MaxNumMultipliers = 4#10
    MaxPos = 4# 10

    dim = rand(1:MaxDim)

    if dim == 2
        mytol = 0.2
    elseif dim==3
        mytol = 0.13
    else
        mytol=0.1
    end

    A = round.(rand(dim, dim), digits = maxdigits)
    while abs(det(A)) < mytol
        A = round.(rand(dim, dim), digits = maxdigits)
    end


    L = alfa.Lattice(A)
    domain =
        [L.dA * unique(round.(rand(dim), digits = maxdigits)) for _ = 1:rand(1:MaxNumElements)]
    if single_domain
        codomain = domain
    else
        codomain = [
            L.dA * round.(rand(dim), digits = maxdigits)
            for _ = 1:rand(1:MaxNumElements)
        ]
    end
    C = alfa.Crystal(L, domain, codomain)
    S = alfa.CrystalOperator(C)

    M = rand(1:MaxNumMultipliers)


    for m in M
        pos = rand(-MaxPos:MaxPos, dim)
        m =
            round.(
                rand(ComplexF64, C.size_codomain, C.size_domain),
                digits = maxdigits,
            )
        push!(S, alfa.Multiplier(pos, m))
    end

    return S
end

function Base.rand(
    A::alfa.CrystalOperator;
    domain_eq_Adomain = false,
    domain_eq_Acodomain = false,
    codomain_eq_Adomain = false,
    codomain_eq_Acodomain = false,
    maxdigits = 3
)
    if A.C.dim > 2
        MaxLatticeSize = 1
    else
        MaxLatticeSize = 2
    end
    MaxNumElements = 2 #10
    MaxNumMultipliers = 4 #10
    MaxPos = 4# 10

    B = alfa.lll(rand(-MaxLatticeSize:MaxLatticeSize, A.C.dim, A.C.dim))
    while abs(det(B)) < 1e-1
        B = alfa.lll(rand(-MaxLatticeSize:MaxLatticeSize, A.C.dim, A.C.dim))
    end
    B = A.C.L.A * B # make it a sublattice of A

    L = alfa.Lattice(B)
    if domain_eq_Adomain || domain_eq_Acodomain || codomain_eq_Acodomain || codomain_eq_Adomain
        C = alfa.wrtLattice(A.C, L)
    end

    if domain_eq_Adomain
        domain = C.Domain
    elseif domain_eq_Acodomain
        domain = C.Codomain
    else
        domain = [
            L.dA * unique(round.(rand(L.dim), digits = maxdigits))
            for _ = 1:rand(1:MaxNumElements)
        ]
    end
    if codomain_eq_Acodomain
        codomain = C.Codomain
    elseif codomain_eq_Adomain
        codomain = C.Domain
    else
        codomain = [
            L.dA * unique(round.(rand(L.dim), digits = maxdigits))
            for _ = 1:rand(1:MaxNumElements)
        ]
    end

    C = alfa.Crystal(L, domain, codomain)
    S = alfa.CrystalOperator(C)

    M = rand(1:MaxNumMultipliers)


    for m in M
        pos = rand(-MaxPos:MaxPos, A.C.dim)
        m =
            round.(
                rand(ComplexF64, C.size_codomain, C.size_domain),
                digits = maxdigits,
            )
        push!(S, alfa.Multiplier(pos, m))
    end

    return S
end

end
