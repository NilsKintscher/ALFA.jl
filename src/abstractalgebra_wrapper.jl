## AbstractAlgebra conversion.
function Base.convert(
    ::Type{Matrix{Int}},
    x::AbstractAlgebra.Generic.Mat{BigInt},
)
    m, n = size(x)
    mat = Int[x[i, j] for i = 1:m, j = 1:n]
    return convert(Matrix{Int}, mat)
end

function Base.convert(
    ::Type{MMatrix{M,N,Int}},
    x::AbstractAlgebra.Generic.Mat{BigInt},
) where N where M
    mat = MMatrix{N,M,Int}([x[i, j] for i=1:N, j=1:M])
    return mat
end

function Base.convert(::Type{AbstractAlgebra.Generic.Mat{BigInt}}, mat)
    ms = AbstractAlgebra.MatrixSpace(AbstractAlgebra.ZZ, size(mat)...)
    mat = ms(convert(Matrix{Int}, round.(mat)))
    return mat
end

#Smith normal form wrapper
function snf_with_transform(mat::Matrix) #{T} where T<:Real)
    mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
    (S, U, V) = AbstractAlgebra.snf_with_transform(mat) # U*mat*V = S
    U = convert(MMatrix{size(U)...,Int}, U)
    V = convert(MMatrix{size(V)...,Int}, V)
    S = convert(MMatrix{size(S)...,Int}, S)
    return S, U, V
end

#Smith normal form wrapper
function snf_with_transform(mat::MMatrix{M,N}) where M where N #{T} where T<:Real)
    mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
    (S, U, V) = AbstractAlgebra.snf_with_transform(mat) # U*mat*V = S
    U = convert(MMatrix{size(U)...,Int}, U)
    V = convert(MMatrix{size(V)...,Int}, V)
    S = convert(MMatrix{size(S)...,Int}, S)
    return S, U, V
end

#Hermite normal form wrapper
function hnf(mat::Matrix)
    # computes the H = mat*U, s.t. H is in HNF and U is unimodular.
    mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
    H = AbstractAlgebra.hnf(transpose(mat)) # H = upper triangular
    H = convert(Matrix{Int}, transpose(H))
    return H
end

function hnf(mat::MMatrix{M,N}) where M where N
    # computes the H = mat*U, s.t. H is in HNF and U is unimodular.
    mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
    H = transpose(AbstractAlgebra.hnf(transpose(mat)))
    H = convert(MMatrix{size(H)...,Int}, H)
    return H
end
