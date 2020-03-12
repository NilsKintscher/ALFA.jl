## AbstractAlgebra conversion.
# function Base.convert(
#     ::Type{Matrix{Int}},
#     x::AbstractAlgebra.Generic.Mat{BigInt},
# )
#     m, n = size(x)
#     mat = Int[x[i, j] for i = 1:m, j = 1:n]
#     return convert(Matrix{Int}, mat)
# end
#
# function Base.convert(
#     ::Type{MMatrix{M,N,Int}},
#     x::AbstractAlgebra.Generic.Mat{BigInt},
# ) where N where M
#     mat = MMatrix{N,M,Int}([x[i, j] for i=1:N, j=1:M])
#     return mat
# end
#
# function Base.convert(::Type{AbstractAlgebra.Generic.Mat{BigInt}}, mat)
#     ms = AbstractAlgebra.MatrixSpace(AbstractAlgebra.ZZ, size(mat)...)
#     mat = ms(convert(Matrix{Int}, round.(mat)))
#     return mat
# end

function Base.convert(::Type{Nemo.fmpz_mat}, mat::Union{Matrix{T}, MMatrix{N,N,T}}) where {N, T<:Real}
    return  Nemo.MatrixSpace(Nemo.ZZ, size(mat)...)( Int.(round.(mat)) )
end

function Base.convert(::Type{MMatrix{N,N,Int}}, mat::Nemo.fmpz_mat) where {N}
    A = Matrix{Int}(mat)
    return convert(MMatrix{N,N,Int}, A)
end

function Base.convert(
    ::Type{Matrix{Int}},
    mat::Nemo.fmpz_mat,
)
    return Matrix{Int}(mat)
end


#Smith normal form wrapper
# function snf_with_transform(mat::Matrix) #{T} where T<:Real)
#     mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
#     (S, U, V) = AbstractAlgebra.snf_with_transform(mat) # U*mat*V = S
#     U = convert(MMatrix{size(U)...,Int}, U)
#     V = convert(MMatrix{size(V)...,Int}, V)
#     S = convert(MMatrix{size(S)...,Int}, S)
#     return S, U, V
# end

#Smith normal form wrapper
# function snf_with_transform2(mat::MMatrix{M,N}) where M where N #{T} where T<:Real)
#     mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
#     (S, U, V) = AbstractAlgebra.snf_with_transform(mat) # U*mat*V = S
#     U = convert(MMatrix{size(U)...,Int}, U)
#     V = convert(MMatrix{size(V)...,Int}, V)
#     S = convert(MMatrix{size(S)...,Int}, S)
#     return S, U, V
# end

function snf_with_transform(mat::MMatrix{M,N}) where M where N #{T} where T<:Real)
    mat = convert(Nemo.fmpz_mat, mat)
    (S, U, V) = Nemo.snf_with_transform(mat) # U*mat*V = S
    U = convert(MMatrix{size(U)...,Int}, U)
    V = convert(MMatrix{size(V)...,Int}, V)
    S = convert(MMatrix{size(S)...,Int}, S)
    return S, U, V
end

function snf_with_transform(mat::Matrix) #where M where N #{T} where T<:Real)
    mat = convert(Nemo.fmpz_mat, mat)
    (S, U, V) = Nemo.snf_with_transform(mat) # U*mat*V = S
    U = convert(Matrix{Int}, U)
    V = convert(Matrix{Int}, V)
    S = convert(Matrix{Int}, S)
    return S, U, V
end

#Hermite normal form wrapper
# function hnf(mat::Matrix)
#     # computes the H = mat*U, s.t. H is in HNF and U is unimodular.
#     mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
#     return mat
#     H = AbstractAlgebra.hnf(transpose(mat)) # H = upper triangular
#     H = convert(Matrix{Int}, transpose(H))
#     return H
# end

# function hnf2(mat::MMatrix{M,N}) where M where N
#     # computes the H = mat*U, s.t. H is in HNF and U is unimodular.
#     mat = convert(AbstractAlgebra.Generic.Mat{BigInt}, mat)
#     H = transpose(AbstractAlgebra.hnf(transpose(mat)))
#     H = convert(MMatrix{size(H)...,Int}, H)
#     return H
# end


function hnf(mat::MMatrix{M,N}) where M where N
    # computes the H = mat*U, s.t. H is in HNF and U is unimodular.
    mat = convert(Nemo.fmpz_mat, mat)
    H = transpose(Nemo.hnf(transpose(mat)))
    H = convert(MMatrix{size(H)...,Int}, H)
    return H
end

function hnf(mat::Matrix)
    # computes the H = mat*U, s.t. H is in HNF and U is unimodular.
    mat = convert(Nemo.fmpz_mat, mat)
    H = transpose(Nemo.hnf(transpose(mat)))
    H = convert(Matrix{Int}, H)
    return H
end

function lll(mat::MMatrix{M,N}) where M where N
    # computes L, such that mat*T = L for some unimodular T
    mat = convert(Nemo.fmpz_mat, mat)
    L = Nemo.lll(transpose(mat))
    L = transpose(L)
    L = convert(MMatrix{size(L)...,Int}, L)
    return L
end

function lll(mat::Matrix)
    # computes L, such that mat*T = L for some unimodular T
    mat = convert(Nemo.fmpz_mat, mat)
    L = Nemo.lll(transpose(mat))
    L = transpose(L)
    L = convert(Matrix{Int}, L)
    return L
end

#
# function lll_with_transform(mat::MMatrix{M,N}) where M where N
#     mat = convert(Nemo.fmpz_mat, mat)
#     (L,T) = Nemo.lll_with_transform(transpose(mat))
#     L = transpose(L)
#     T = transpose(T)
#     L = convert(MMatrix{size(L)...,Int}, L)
#     T = convert(MMatrix{size(T)...,Int}, T)
#     return L, T
# end
