#
module gallery
using ..alfa
using LinearAlgebra

function Laplace(;N = 2, h=1, T = Float64)
    if N == 1
        A = h*[1]
        Domain = [0]
        Codomain = [0]

        C = alfa.Crystal{N,T}(A, Domain, Codomain)
        L = alfa.CrystalOperator{N,T}(C)

        push!(L, alfa.Multiplier([0], [-2/h^2]))
        push!(L, alfa.Multiplier([-1], [1/h^2]))
        push!(L, alfa.Multiplier([1], [1/h^2]))
        return L

    elseif N == 2
        A = h*[1 0; 0 1]
        Domain = [0 0]
        Codomain = [0 0]

        C = alfa.Crystal{N,T}(A, Domain, Codomain)
        L = alfa.CrystalOperator{N,T}(C)

        push!(L, alfa.Multiplier([0 0], [-4/h^2]))
        push!(L, alfa.Multiplier([0 -1], [1/h^2]))
        push!(L, alfa.Multiplier([0 1], [1/h^2]))
        push!(L, alfa.Multiplier([1 0], [1/h^2]))
        push!(L, alfa.Multiplier([-1 0], [1/h^2]))
        return L
    end
end

function fw_restriction(;m=1, N = 2, T = Float64)
    if N==1
        A = 2^m * [1]
        Domain = [[0], [1]]
        Codomain = [0]

        C = alfa.Crystal{1,T}(A, Domain, Codomain)
        L = alfa.CrystalOperator{1,T}(C)

        push!(L, alfa.Multiplier([0], [1 1//2]))
        push!(L, alfa.Multiplier([-1], [0 1//2]))
        return L
    elseif N == 2
        A = 2^m * [1 0; 0 1]
        Domain = [[0, 0], [0, 1], [1, 0], [1, 1]]
        Codomain = [0, 0]

        C = alfa.Crystal{2,T}(A, Domain, Codomain)
        L = alfa.CrystalOperator{2,T}(C)

        push!(L, alfa.Multiplier([0 0], [1 1 // 2 1 // 2 1 // 4]))
        push!(L, alfa.Multiplier([-1 0], [0 0 1 // 2 1 // 4]))
        push!(L, alfa.Multiplier([0 -1], [0 1 // 2 0 1 // 4]))
        push!(L, alfa.Multiplier([-1 -1], [0 0 0 1 // 4]))
        return L
    end
end

function graphene_tight_binding(;t = nothing, T=Float64)
    if t == nothing
        t = [0, -1, 0, 0]
    end
    A = 1 // 2 * [3 3; sqrt(3) -sqrt(3)]
    Domain = [j // 3 * A * [1, 1] for j in [1, 2]]

    C = alfa.Crystal{2,T}(A, Domain)
    L = alfa.CrystalOperator{2,T}(C)

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

function graphene_dirac_restriction(;T=Float64, m = 1, wl = 0.5, wlh = 0)
    ws = wl + wlh - 1

    A = 2^(m - 1) * 1 // 2 * [3 3; sqrt(3) -sqrt(3)]
    Ac = 2 * A

    s01 = [j // 3 * A * [1, 1] for j in [1, 2]]

    Domain = [A * x + y for x in [[0, 0], [1, 0], [0, 1], [1, 1]] for y in s01]
    Codomain = [j // 3 * Ac * [1, 1] for j in [1, 2]]

    C = alfa.Crystal{2,T}(Ac, Domain, Codomain)
    L = alfa.CrystalOperator{2,T}(C)

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



function Base.rand(
    ::Type{alfa.CrystalOperator{N,T}};
    single_domain = false,
) where {N,T}
    maxdigits = 2

    MaxNumElements = 2 #10
    MaxNumMultipliers = 4#10
    MaxPos = 4# 10

    if N == 2
        mytol = 0.2
    elseif N == 3
        mytol = 0.13
    else
        mytol = 0.1
    end

    num = rand(1:10^maxdigits, N, N)
    den = rand(1:10^maxdigits, N, N)
    if T <: Rational
        X = T.parameters[1]
        A = [X(x) // X(y) for (x, y) in zip(num, den)]
        while abs(det(A)) < mytol
            num = rand(1:10^maxdigits, N, N)
            den = rand(1:10^maxdigits, N, N)
            A = [X(x) // X(y) for (x, y) in zip(num, den)]
        end
    else
        A = T.(round.(rand(N, N), digits = maxdigits))
        while abs(det(A)) < mytol
            A = T.(round.(rand(N, N), digits = maxdigits))
        end
    end


    L = alfa.Lattice{N,T}(A)
    if T <: Rational
        denom = [
            [rand(1:10^maxdigits) for _ = 1:N] for _ = 1:rand(1:MaxNumElements)
        ]
        myrand(x) = rand(1:x) // x
        domain = unique([
            L.A * [myrand(rand(1:10^maxdigits)) for _ = 1:N] for _ = 1:rand(1:MaxNumElements)
        ])
    else
        domain = unique([
            L.A * round.(rand(N), digits = maxdigits)
            for _ = 1:rand(1:MaxNumElements)
        ])
    end
    if single_domain
        codomain = domain
    else
        if T <: Rational
            codomain = unique([
                L.A * [myrand(rand(1:10^maxdigits)) for _ = 1:N] for _ = 1:rand(1:MaxNumElements)
            ])
        else
            codomain = unique([
                L.A * round.(rand(N), digits = maxdigits)
                for _ = 1:rand(1:MaxNumElements)
            ])
        end
    end
    C = alfa.Crystal{N,T}(L, domain, codomain)
    S = alfa.CrystalOperator{N,T}(C)

    M = rand(1:MaxNumMultipliers)


    for m in 1:M
        pos = rand(-MaxPos:MaxPos, N)
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
    A::alfa.CrystalOperator{N,T};
    domain_eq_Adomain = false,
    domain_eq_Acodomain = false,
    codomain_eq_Adomain = false,
    codomain_eq_Acodomain = false,
    maxdigits = 2,
) where {N,T}

    MaxLatticeSize = 2
    MaxNumElements = 2 #10
    MaxNumMultipliers = 4 #10
    MaxPos = 4# 10

    B = alfa.lll(rand(-MaxLatticeSize:MaxLatticeSize, N, N))
    while abs(det(B)) > 8 || abs(det(B)) < 1
        B = alfa.lll(rand(-MaxLatticeSize:MaxLatticeSize, N, N))
    end
    B = A.C.L.A * B # make it a sublattice of A

    L = alfa.Lattice{N,T}(B)
    if domain_eq_Adomain ||
       domain_eq_Acodomain || codomain_eq_Acodomain || codomain_eq_Adomain
        C = alfa.wrtLattice(A.C, L)
    end

    myrand(x) = rand(1:x) // x

    if domain_eq_Adomain
        domain = C.Domain
    elseif domain_eq_Acodomain
        domain = C.Codomain
    else

        if T <: Rational
            domain = unique([
                L.A * [myrand(rand(1:10^maxdigits)) for _ = 1:L.dim] for _ = 1:rand(1:MaxNumElements)
            ])
        else
            domain = unique([
                L.A * round.(rand(N), digits = maxdigits)
                for _ = 1:rand(1:MaxNumElements)
            ])
        end
    end
    if codomain_eq_Acodomain
        codomain = C.Codomain
    elseif codomain_eq_Adomain
        codomain = C.Domain
    else
        if T <: Rational
            codomain = unique([
                L.A * [myrand(rand(1:10^maxdigits)) for _ = 1:L.dim] for _ = 1:rand(1:MaxNumElements)
            ])
        else
            codomain = unique([
                L.A * round.(rand(N), digits = maxdigits)
                for _ = 1:rand(1:MaxNumElements)
            ])
        end
    end

    C = alfa.Crystal{N,T}(L, domain, codomain)
    S = alfa.CrystalOperator{N,T}(C)

    M = rand(1:MaxNumMultipliers)


    for m in 1:M
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
