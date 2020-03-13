
function Base.convert(::Type{Nemo.fmpz_mat}, mat::Union{Matrix{T}, MMatrix{N,N,T}}) where {N, T<:Real}
    return  Nemo.MatrixSpace(Nemo.ZZ, size(mat)...)( BigInt.(round.(mat)) )
end

function Base.convert(::Type{MMatrix{N,N,BigInt}}, mat::Nemo.fmpz_mat) where {N}
    A = Matrix{BigInt}(mat)
    return convert(MMatrix{N,N,BigInt}, A)
end

function Base.convert(
    ::Type{Matrix{BigInt}},
    mat::Nemo.fmpz_mat,
)
    return Matrix{BigInt}(mat)
end



function snf_with_transform(mat::MMatrix{M,N}) where M where N #{T} where T<:Real)
    mat = convert(Nemo.fmpz_mat, mat)
    (S, U, V) = Nemo.snf_with_transform(mat) # U*mat*V = S
    U = convert(MMatrix{size(U)...,BigInt}, U)
    V = convert(MMatrix{size(V)...,BigInt}, V)
    S = convert(MMatrix{size(S)...,BigInt}, S)
    return S, U, V
end

function snf_with_transform(mat::Matrix) #where M where N #{T} where T<:Real)
    mat = convert(Nemo.fmpz_mat, mat)
    (S, U, V) = Nemo.snf_with_transform(mat) # U*mat*V = S
    U = convert(Matrix{BigInt}, U)
    V = convert(Matrix{BigInt}, V)
    S = convert(Matrix{BigInt}, S)
    return S, U, V
end


function hnf(mat::MMatrix{M,N}) where M where N
    # computes the H = mat*U, s.t. H is in HNF and U is unimodular.
    mat = convert(Nemo.fmpz_mat, mat)
    H = transpose(Nemo.hnf(transpose(mat)))
    H = convert(MMatrix{size(H)...,BigInt}, H)
    return H
end

function hnf(mat::Matrix)
    # computes the H = mat*U, s.t. H is in HNF and U is unimodular.
    mat = convert(Nemo.fmpz_mat, mat)
    H = transpose(Nemo.hnf(transpose(mat)))
    H = convert(Matrix{BigInt}, H)
    return H
end

function lll(mat::MMatrix{M,N}) where M where N
    # computes L, such that mat*T = L for some unimodular T
    mat = convert(Nemo.fmpz_mat, mat)
    L = Nemo.lll(transpose(mat))
    L = transpose(L)
    L = convert(MMatrix{size(L)...,BigInt}, L)
    return L
end

function lll(mat::Matrix)
    # computes L, such that mat*T = L for some unimodular T
    mat = convert(Nemo.fmpz_mat, mat)
    L = Nemo.lll(transpose(mat))
    L = transpose(L)
    L = convert(Matrix{BigInt}, L)
    return L
end
