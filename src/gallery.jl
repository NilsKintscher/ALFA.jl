#
module gallery
using ..alfa
using LinearAlgebra

function Laplace(N = 2, T = Float64)
    if N == 1
        A = [1]
        Domain = [0]
        Codomain = [0]

        C = alfa.Crystal{N,T}(A, Domain, Codomain)
        L = alfa.CrystalOperator{N,T}(C)

        push!(L, alfa.Multiplier([0], [-2]))
        push!(L, alfa.Multiplier([-1], [1]))
        push!(L, alfa.Multiplier([1], [1]))
        return L

    elseif N == 2
        A = [1 0; 0 1]
        Domain = [0 0]
        Codomain = [0 0]

        C = alfa.Crystal{N,T}(A, Domain, Codomain)
        L = alfa.CrystalOperator{N,T}(C)

        push!(L, alfa.Multiplier([0 0], [-4]))
        push!(L, alfa.Multiplier([0 -1], [1]))
        push!(L, alfa.Multiplier([0 1], [1]))
        push!(L, alfa.Multiplier([1 0], [1]))
        push!(L, alfa.Multiplier([-1 0], [1]))
        return L
    end
end

function fw_restriction(N = 2, T = Float64)
    if N == 2
        A = 2 * [1 0; 0 1]
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



function Base.rand(
    ::Type{alfa.CrystalOperator{N,T}};
    single_domain = false,
) where {N,T}
    maxdigits = 2

    MaxNumElements = 2 #10
    MaxNumMultipliers = 1#10
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


    for m in M
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
    MaxNumMultipliers = 1 #10
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












#
#
#
# A= alfa.CrystalOperator{2,Float64}(alfa.Crystal{2,Float64}(alfa.Lattice{2,Float64}([-0.52 -0.7; -0.32 -0.85]), SArray{Tuple{2},Float64,1,2}[[1.0408, 1.0178], [0.54, 0.391]], SArray{Tuple{2},Float64,1,2}[[0.511, 0.40249999999999997], [0.7726, 0.6641]], false), SortedSet(alfa.Multiplier[alfa.Multiplier{2}([2, 2], Complex{Float64}[0.74 + 0.09im 0.93 + 0.76im; 0.21 +
# 0.79im 0.61 + 0.54im])],Base.Order.ForwardOrdering()), false)
#
#  B= alfa.CrystalOperator{2,Float64}(alfa.Crystal{2,Float64}(alfa.Lattice{2,Float64}([1.4 -1.74; 1.7 -1.49]), SArray{Tuple{2},Float64,1,2}[[1.0408, 1.0178], [0.54, 0.391], [1.5608, 1.3378], [1.06, 0.7110000000000001], [1.7408, 1.8678], [1.24, 1.241], [2.2607999999999997, 2.1878], [1.76, 1.561]], SArray{Tuple{2},Float64,1,2}[[0.511, 0.40249999999999997], [0.7726, 0.6641], [1.0310000000000001, 0.7224999999999999], [1.2926, 0.9841], [1.2109999999999999, 1.2525], [1.4726, 1.5141], [1.7309999999999999, 1.5724999999999998], [1.9926, 1.8340999999999998]], false), SortedSet(alfa.Multiplier[alfa.Multiplier{2}([-1, 0], Complex{Float64}[0.64 + 0.4im 0.93 + 0.77im 0.38 + 0.02im 0.15 + 0.54im 0.41 + 0.19im 0.67 + 0.56im 0.15 + 0.31im 0.01 + 0.05im; 0.41 + 0.57im 0.85 + 0.76im 0.98 + 0.42im 0.55 + 0.09im 0.96 + 1.0im 0.51 + 0.53im 0.35 + 0.68im 0.07 + 0.6im; 0.51 + 0.7im 0.3 + 0.85im 0.46 + 0.48im 0.89 + 0.45im 0.65 + 0.36im 0.36 + 0.94im 0.84 + 0.35im 0.38 + 0.92im; 0.57 + 0.77im 0.09 + 0.34im 0.25 + 0.96im 0.44 + 0.56im 0.11 + 0.45im 0.28 + 0.59im 0.79 + 0.22im 0.16 + 0.96im; 0.18 + 0.75im 0.36 + 0.18im 0.98 + 0.32im 0.41 + 0.6im 0.39 + 0.5im 0.25 + 0.19im 0.56 + 0.15im 0.81 + 0.4im; 0.73 + 0.8im 0.07 + 0.19im 0.42 + 0.28im 0.96 + 0.53im 0.56 + 0.77im 0.04 + 0.92im 0.3 + 0.09im 0.09 + 0.31im; 0.01 + 0.07im 0.57 + 0.91im 0.39 + 0.25im 0.37 + 0.66im 0.0 + 0.48im 0.85 + 0.5im 0.64 + 0.87im 0.35 +
# 0.89im; 0.2 + 0.93im 0.16 + 0.48im 0.62 + 0.03im 0.91 + 0.84im 0.62 + 0.49im 0.42 + 0.14im 0.77 + 0.26im 0.93 + 0.28im])],
# Base.Order.ForwardOrdering()), false)
#
# AB1 = A+B
#
# BA1 = B+A
#
# (A1, B1) = alfa.wrtSameLatticeAndNormalize(AB1,BA1)
#
# M1 = collect(A1.M)
# M2 = collect(B1.M)



O= alfa.CrystalOperator{2,Float64}(alfa.Crystal{2,Float64}(alfa.Lattice{2,Float64}([0.45 0.95; 0.62 0.28]), SArray{Tuple{2},Float64,1,2}[[0.4675, 0.37660000000000005]], SArray{Tuple{2},Float64,1,2}[[1.2315, 0.8016000000000001], [1.201, 0.6464000000000001]], false), SortedSet(alfa.Multiplier[alfa.Multiplier{2}([0, 0], Complex{Float64}[0.86 + 0.92im; 0.32 + 0.35im])],
Base.Order.ForwardOrdering()), false)
C=alfa.CrystalOperator{2,Float64}(alfa.Crystal{2,Float64}(alfa.Lattice{2,Float64}([0.9 -2.35; 1.24 -1.1800000000000002]), SArray{Tuple{2},Float64,1,2}[[0.4675, 0.37660000000000005], [1.4175, 0.6566000000000001], [2.3674999999999997, 0.9366000000000001], [3.3175, 1.2166000000000001]], SArray{Tuple{2},Float64,1,2}[[1.2315, 0.8016000000000001], [1.201, 0.6464000000000001], [2.1814999999999998, 1.0816000000000001], [2.151, 0.9264000000000001], [3.1315, 1.3616000000000001], [3.101, 1.2064000000000001], [4.0815, 1.6416000000000002], [4.051, 1.4864000000000002]], false), SortedSet(alfa.Multiplier[alfa.Multiplier{2}([-3, 0], Complex{Float64}[1.0 + 0.92im 0.86 + 0.78im 0.09 + 0.07im 1.0 + 0.8im; 0.36 + 0.71im 0.27 + 0.19im 0.02 + 0.87im 0.27 + 1.0im; 0.43 + 0.33im 0.13 + 0.02im 0.41 + 0.62im 0.51 + 0.44im; 0.12 + 0.25im 0.84 + 0.57im 0.2 + 0.26im 0.38 + 0.14im; 0.96 + 0.49im 0.3 + 0.62im 0.64 + 0.47im 0.43 + 0.96im; 0.24 + 0.0im 0.06 + 0.09im 0.7 + 0.01im 0.22 + 0.83im; 0.84 + 0.95im 0.48 + 0.11im 0.32 + 0.52im 0.12 + 0.7im; 0.21 + 1.0im 0.22 + 0.78im 0.86 + 0.66im 0.9 + 0.08im])],
Base.Order.ForwardOrdering()), false)

A=alfa.CrystalOperator{2,Float64}(alfa.Crystal{2,Float64}(alfa.Lattice{2,Float64}([2.35 0.04999999999999993; 1.1800000000000002 -0.96]), SArray{Tuple{2},Float64,1,2}[[1.2315, 0.8016000000000001], [1.201, 0.6464000000000001], [2.1814999999999998, 1.0816000000000001], [2.151, 0.9264000000000001], [3.1315, 1.3616000000000001], [3.101, 1.2064000000000001], [4.0815, 1.6416000000000002], [4.051, 1.4864000000000002], [5.031499999999999, 1.9216000000000002], [5.0009999999999994, 1.7664000000000002]], SArray{Tuple{2},Float64,1,2}[[2.35, 0.7170000000000002]], false), SortedSet(alfa.Multiplier[alfa.Multiplier{2}([2, -3],
Complex{Float64}[0.52 + 0.95im 0.99 + 0.21im 0.77 + 0.46im 0.46 + 0.76im 0.73 + 0.85im 0.83 + 0.78im 0.74 + 0.27im 0.17 + 0.86im 0.66 + 0.44im 0.86 + 0.44im])],
Base.Order.ForwardOrdering()), false)
