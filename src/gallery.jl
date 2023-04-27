#
module gallery
using ..ALFA
using LinearAlgebra

"""
    Laplace(;N = 2, h=1, T = Float64)

Creates a CrystalOperator corresponding to the (central differences) discretization of the Laplace operator ``Δ=\\sum_j \\frac{\\partial^2}{\\partial x_j^2} ``in N dimensions.
- h corresponds to the distance of the grid points.

"""
function Laplace(;N = 2, h=1, T = Float64)
    if N == 1
        A = h*[1]
        Domain = [0]
        Codomain = [0]

        C = ALFA.Crystal{N,T}(A, Domain, Codomain)
        L = ALFA.CrystalOperator{N,T}(C)

        push!(L, ALFA.Multiplier([0], [-2/h^2]))
        push!(L, ALFA.Multiplier([-1], [1/h^2]))
        push!(L, ALFA.Multiplier([1], [1/h^2]))
        return L

    elseif N == 2
        A = h*[1 0; 0 1]
        Domain = [0 0]
        Codomain = [0 0]

        C = ALFA.Crystal{N,T}(A, Domain, Codomain)
        L = ALFA.CrystalOperator{N,T}(C)

        push!(L, ALFA.Multiplier([0 0], [-4/h^2]))
        push!(L, ALFA.Multiplier([0 -1], [1/h^2]))
        push!(L, ALFA.Multiplier([0 1], [1/h^2]))
        push!(L, ALFA.Multiplier([1 0], [1/h^2]))
        push!(L, ALFA.Multiplier([-1 0], [1/h^2]))
        return L
    end
end

"""
    fw_restriction(;m=1, N = 2, T = Float64)

Full weighting restriction operator in N dimensions.

"""
function fw_restriction(;m=1, N = 2, T = Float64)
    if N==1
        A = 2^m * [1]
        Domain = [[0], [1]]
        Codomain = [0]

        C = ALFA.Crystal{1,T}(A, Domain, Codomain)
        L = ALFA.CrystalOperator{1,T}(C)

        push!(L, ALFA.Multiplier([0], [1 1//2]))
        push!(L, ALFA.Multiplier([-1], [0 1//2]))
        return L
    elseif N == 2
        A = 2^m * [1 0; 0 1]
        Domain = [[0, 0], [0, 1], [1, 0], [1, 1]]
        Codomain = [0, 0]

        C = ALFA.Crystal{2,T}(A, Domain, Codomain)
        L = ALFA.CrystalOperator{2,T}(C)

        push!(L, ALFA.Multiplier([0 0], [1 1 // 2 1 // 2 1 // 4]))
        push!(L, ALFA.Multiplier([-1 0], [0 0 1 // 2 1 // 4]))
        push!(L, ALFA.Multiplier([0 -1], [0 1 // 2 0 1 // 4]))
        push!(L, ALFA.Multiplier([-1 -1], [0 0 0 1 // 4]))
        return L
    end
end

"""
    graphene_tight_binding(;t = nothing, T=Float64)

Creates a CrystalOperator corresponding to the tight-binding Hamiltonian of graphene.
The vector t should contain the hopping parameter.
See [Section 6.1, 1] Kahl, K., Kintscher, N. Automated local Fourier analysis (ALFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>

"""
function graphene_tight_binding(;t = nothing, T=Float64)
    if t === nothing
        t = [0, -1, 0, 0]
    end
    A = 1 // 2 * [3 3; sqrt(3) -sqrt(3)]
    Domain = [j // 3 * A * [1, 1] for j in [1, 2]]

    C = ALFA.Crystal{2,T}(A, Domain)
    L = ALFA.CrystalOperator{2,T}(C)

    push!(L, ALFA.Multiplier([0, 0], [t[1] t[2]; t[2] t[1]]))

    push!(L, ALFA.Multiplier([-1, 0], [t[3] t[2]; 0 t[3]]))
    push!(L, ALFA.Multiplier([0, -1], [t[3] t[2]; 0 t[3]]))

    push!(L, ALFA.Multiplier([1, 0], [t[3] 0; t[2] t[3]]))
    push!(L, ALFA.Multiplier([0, 1], [t[3] 0; t[2] t[3]]))

    push!(L, ALFA.Multiplier([1, -1], [t[3] t[4]; t[4] t[3]]))
    push!(L, ALFA.Multiplier([-1, 1], [t[3] t[4]; t[4] t[3]]))

    push!(L, ALFA.Multiplier([1, 1], [0 0; t[4] 0]))
    push!(L, ALFA.Multiplier([-1, -1], [0 t[4]; 0 0]))

    return L
end

"""
    graphene_dirac_restriction(;t = nothing, T=Float64)

Creates a restriction CrystalOperator that conserveres the dirac function.
See [Section 6.1, 1] Kahl, K., Kintscher, N. Automated local Fourier analysis (ALFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>

"""
function graphene_dirac_restriction(;T=Float64, m = 1, wl = 0.5, wlh = 0)
    ws = wl + wlh - 1

    A = 2^(m - 1) * 1 // 2 * [3 3; sqrt(3) -sqrt(3)]
    Ac = 2 * A

    s01 = [j // 3 * A * [1, 1] for j in [1, 2]]

    Domain = [A * x + y for x in [[0, 0], [1, 0], [0, 1], [1, 1]] for y in s01]
    Codomain = [j // 3 * Ac * [1, 1] for j in [1, 2]]

    C = ALFA.Crystal{2,T}(Ac, Domain, Codomain)
    L = ALFA.CrystalOperator{2,T}(C)

    m = zeros(2, 8)
    m[2, 1] = wl
    push!(L, ALFA.Multiplier([1, 1], m))


    m = zeros(2, 8)
    m[2, [1, 3, 5]] = [ws, wlh, ws]
    push!(L, ALFA.Multiplier([1, 0], m))

    m = zeros(2, 8)
    m[2, [1, 3, 5]] = [ws, ws, wlh]
    push!(L, ALFA.Multiplier([0, 1], m))

    m = zeros(2, 8)
    m[2, 5] = wl
    m[1, 6] = wl
    push!(L, ALFA.Multiplier([1, -1], m))

    m = zeros(2, 8)
    m[1, [2, 4, 6, 8]] = [1, ws, ws, wlh]
    m[2, [1, 3, 5, 7]] = [wlh, ws, ws, 1]
    push!(L, ALFA.Multiplier([0, 0], m))

    m = zeros(2, 8)
    m[2, 3] = wl
    m[1, 4] = wl
    push!(L, ALFA.Multiplier([-1, 1], m))

    m = zeros(2, 8)
    m[1, [4, 6, 8]] = [wlh, ws, ws]
    push!(L, ALFA.Multiplier([0, -1], m))

    m = zeros(2, 8)
    m[1, [4, 6, 8]] = [ws, wlh, ws]
    push!(L, ALFA.Multiplier([-1, 0], m))

    m = zeros(2, 8)
    m[1, 8] = wl
    push!(L, ALFA.Multiplier([-1, -1], m))

    return L
end

"""
    curlcurl(sigma; T=Float64)

Creates the curlcurl operator.
See [Section 6.2, 1] Kahl, K., Kintscher, N. Automated local Fourier analysis (ALFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>

"""
function curlcurl(sigma=0; T=Float64)
    A = [1 0; 0 1]
    Domain = [[.5, 0],[0, .5]]

    C = ALFA.Crystal{2,T}(A, Domain)
    L = ALFA.CrystalOperator{2,T}(C)

    push!(L, ALFA.Multiplier([0 0], [2+2/3*sigma -1; -1 2+2/3*sigma]))
    push!(L, ALFA.Multiplier([0 -1], [-1+sigma/6 1; 0 0]))
    push!(L, ALFA.Multiplier([0 1], [-1+sigma/6 0; 1 0]))
    push!(L, ALFA.Multiplier([-1 0], [0 0;1 -1+sigma/6]))
    push!(L, ALFA.Multiplier([1 0], [0 1;0 -1+sigma/6]))
    push!(L, ALFA.Multiplier([-1 1], [0 0;-1 0]))
    push!(L, ALFA.Multiplier([1 -1], [0 -1;0 0]))
    return L
end

"""
    curlcurl_restriction(sigma; T=Float64)

Creates the restriction operator corresponding to the curl curl operator.
See [Section 6.2, 1] Kahl, K., Kintscher, N. Automated local Fourier analysis (ALFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>

"""
function curlcurl_restriction(; T=Float64)
    A = 2*[1 0; 0 1]
    s = [[.5, 0],[0, .5]]
    Domain = [sj+[x,y] for sj in s for x in [0,1] for y in [0,1]]
    Codomain = [2*sj for sj in s]

    C = ALFA.Crystal{2,T}(A, Domain, Codomain)
    L = ALFA.CrystalOperator{2,T}(C)
    L = ALFA.normalize(L)
    #
    m = zeros(2,8)
    m[1, [1,2]] .= .5
    m[1, [5,6]] .= .25
    m[2, [3,7]] .= .5
    m[2, [4,8]] .= .25
    push!(L, ALFA.Multiplier([0 0], m))

    m = zeros(2,8)
    m[2, [4,8]] .= .25
    push!(L, ALFA.Multiplier([0 -1], m))

    m = zeros(2,8)
    m[1, [5,6]] .= .25
    push!(L, ALFA.Multiplier([-1 0], m))

    return L
end



function Base.rand(
    ::Type{ALFA.CrystalOperator{N,T}};
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


    L = ALFA.Lattice{N,T}(A)
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
    C = ALFA.Crystal{N,T}(L, domain, codomain)
    S = ALFA.CrystalOperator{N,T}(C)

    M = rand(1:MaxNumMultipliers)


    for m in 1:M
        pos = rand(-MaxPos:MaxPos, N)
        m =
            round.(
                rand(ComplexF64, C.size_codomain, C.size_domain),
                digits = maxdigits,
            )
        push!(S, ALFA.Multiplier(pos, m))
    end

    return S
end

function Base.rand(
    A::ALFA.CrystalOperator{N,T};
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

    B = ALFA.lll(rand(-MaxLatticeSize:MaxLatticeSize, N, N))
    while abs(det(B)) > 8 || abs(det(B)) < 1
        B = ALFA.lll(rand(-MaxLatticeSize:MaxLatticeSize, N, N))
    end
    B = A.C.L.A * B # make it a sublattice of A

    L = ALFA.Lattice{N,T}(B)
    if domain_eq_Adomain ||
       domain_eq_Acodomain || codomain_eq_Acodomain || codomain_eq_Adomain
        C = ALFA.wrtLattice(A.C, L)
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

    C = ALFA.Crystal{N,T}(L, domain, codomain)
    S = ALFA.CrystalOperator{N,T}(C)

    M = rand(1:MaxNumMultipliers)


    for m in 1:M
        pos = rand(-MaxPos:MaxPos, A.C.dim)
        m =
            round.(
                rand(ComplexF64, C.size_codomain, C.size_domain),
                digits = maxdigits,
            )
        push!(S, ALFA.Multiplier(pos, m))
    end

    return S
end

end
