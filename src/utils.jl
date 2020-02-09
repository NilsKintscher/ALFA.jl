## AbstractAlgebra conversion.
function Base.convert(
    ::Type{Matrix{Int}},
    x::AbstractAlgebra.Generic.Mat{BigInt},
)
    m, n = size(x)
    mat = Int[x[i, j] for i = 1:m, j = 1:n]
    return convert(Matrix{Int}, mat)
end

function Base.convert(::Type{AbstractAlgebra.Generic.Mat{BigInt}}, mat::Matrix)
    ms = AbstractAlgebra.MatrixSpace(AbstractAlgebra.ZZ, size(mat)...)
    mat = ms(convert(Matrix{Int}, round.(mat)))
    return mat
end

#Smith normal form wrapper
function snf_with_transform(mat::Matrix) #{T} where T<:Real)
    mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
    (S, U, V) = AbstractAlgebra.snf_with_transform(mat) # U*mat*V = S
    U = convert(Matrix{Int}, U)
    V = convert(Matrix{Int}, V)
    S = convert(Matrix{Int}, S)
    return S, U, V
end

#Hermite normal form wrapper
function hnf(mat::Matrix)
    mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
    H = AbstractAlgebra.hnf(mat) # H = upper triangular
    H = convert(Matrix{Int}, H)
    return H
end
