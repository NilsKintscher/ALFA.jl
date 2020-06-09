
function Base.convert(
    ::Type{Nemo.fmpz_mat},
    mat::Union{Matrix{T},MMatrix{N,N,T}},
) where {N,T<:Real}
    return Nemo.MatrixSpace(Nemo.ZZ, size(mat)...)(BigInt.(round.(mat)))
end

function Base.convert(::Type{MMatrix{N,N,BigInt}}, mat::Nemo.fmpz_mat) where {N}
    A = Matrix{BigInt}(mat)
    return convert(MMatrix{N,N,BigInt}, A)
end

function Base.convert(::Type{Matrix{BigInt}}, mat::Nemo.fmpz_mat)
    return Matrix{BigInt}(mat)
end


"""
    snf_with_transform(L::Lattice)
    snf_with_transform(mat::MMatrix{M,N}) where {M,N}
    snf_with_transform(mat::Matrix)


Wrapper of Nemo.snf_with_transform. Input is converted to BigInt.
Returns (S,U,V) such that U*mat*V = S, where S is the Smith normal form of mat.

# Example
```jldoctest
julia> using ALFA # hide

julia> using LinearAlgebra # hide

julia> mat = rand(1:10, 10, 10);

julia> (S,U,V) = ALFA.snf_with_transform(mat);

julia> norm(U*mat*V - S) ≈ 0
true

julia> abs(det(U)) ≈ abs(det(V)) ≈ 1
true

julia> norm(diagm(diag(S)) - S) ≈ 0
true
```
"""
function snf_with_transform(mat::MMatrix{M,N}) where {M,N}
    mat = convert(Nemo.fmpz_mat, mat)
    (S, U, V) = Nemo.snf_with_transform(mat) # U*mat*V = S
    U = convert(MMatrix{size(U)...,BigInt}, U)
    V = convert(MMatrix{size(V)...,BigInt}, V)
    S = convert(MMatrix{size(S)...,BigInt}, S)
    return S, U, V
end

function snf_with_transform(mat::Matrix)
    mat = convert(Nemo.fmpz_mat, mat)
    (S, U, V) = Nemo.snf_with_transform(mat) # U*mat*V = S
    U = convert(Matrix{BigInt}, U)
    V = convert(Matrix{BigInt}, V)
    S = convert(Matrix{BigInt}, S)
    return S, U, V
end


"""
    hnf(mat::MMatrix{M,N}) where {M, N}
    hnf(mat::Matrix)

Wrapper of Nemo.hnf. Input is converted to BigInt.
Returns H = mat*U, s.t. H is in Hermite Normal Form and U is unimodular.

# Example
```jldoctest
julia> using ALFA # hide

julia> using LinearAlgebra # hide

julia> mat = rand(1:1000, 2, 2);

julia> H = ALFA.hnf(mat);

julia> norm(LinearAlgebra.tril(H) - H) ≈ 0
true
julia> round(abs(det(inv(mat)*H)), digits=5)
1.0
```
"""
function hnf(mat::MMatrix{M,N}) where {M, N}
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

"""
    lll(mat::MMatrix{M,N}) where {M, N}
    lll(mat::Matrix)


Wrapper of Nemo.lll. Input is converted to BigInt.
Applies the LLL-Algorithm to the input mat.
Computes output L, such that mat*T=L for some unimodular T.

# Example
```jldoctest
julia> using ALFA # hide

julia> using LinearAlgebra # hide

julia> mat = rand(1:1000, 2, 2);

julia> L = ALFA.lll(mat);

julia> round(abs(det(inv(mat)*L)), digits=5)
1.0
```
"""
function lll(mat::MMatrix{M,N}) where {M, N}
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
