"""
    symbol(S::CrystalOperator, k; π = π)
    symbol(O::OperatorComposition, k; π = π)

Returns the symbol of the CrystalOperator/OperatorComposition for a given frequency/wavevector k.


# Example
```jldoctest
julia> using ALFA

julia> L = ALFA.gallery.Laplace(N=2);

julia> oc = ALFA.OperatorComposition(:(3*\$L));

julia> ALFA.symbol(oc,[0.5, 0.5]) ≈ [-24] # as f is the identitity, the fourier transform is 1 for all frequencies.
true

```
"""
function symbol(S::CrystalOperator, k; π = π)
    if length(S.M) == 0
        mat = zeros(Int, S.C.size_codomain, S.C.size_domain)
    else
        mat = zeros(eltype(first(S.M).mat), S.C.size_codomain, S.C.size_domain)
        for m in S.M
            mat += m.mat * exp(im * 2π * ((S.C.A * m.pos)'*k))
        end
    end
    return mat
end
symbol(O::OperatorComposition, k; π = π) = O.sym(k,π=π)

"""
    eigvals(
        S::X,
        k::T;
        by = abs,
    ) where {
        T<:Union{AbstractVector,Tuple},
        X<:Union{CrystalOperator,OperatorComposition},
    }

Computes the eigenvalues of the symbol of a CrystalOperator/OperatorComposition wrt frequency k.

The eigenvalues are sorted by `by`.
"""
function eigvals(
    S::X,
    k::T;
    by = abs,
) where {
    T<:Union{AbstractVector,Tuple},
    X<:Union{CrystalOperator,OperatorComposition},
}
    symb = ComplexF64.(symbol(S, k))
    ev = LinearAlgebra.eigvals(symb)
    if by == nothing
        return ev
    else
        return sort(ev, by = by)
    end
end


"""
    eigvals(
        S::X;
        N = 20,
        by = abs,
        unique = false,
        digits = 5,
    ) where {X<:Union{CrystalOperator,OperatorComposition}}

Computes the eigenvalues of the symbol of a CrystalOperator/OperatorComposition.
The Frequency space is divided into N^dim equidistant (unique) points frequencies k.

"""
function eigvals(
    S::X;
    N = 20,
    by = abs,
    unique = false,
    digits = 5,
) where {X<:Union{CrystalOperator,OperatorComposition}}
    k_iter = Iterators.product([
        range(0, stop = 1, length = N + 1)[1:end-1] for _ = 1:S.C.dim
    ]...)
    #kv = reshape([[x...] for x in collect(k_iter)],:)
    #dAkv = [S.C.dA*x for x in kv]
    if unique
        SS = Set{Complex}()
        for k in k_iter
            Λ = Base.unique(round.(
                eigvals(S, S.C.dA * [k...], by = nothing),
                digits = digits,
            ))
            for λ in Λ
                push!(SS, λ)
            end
        end
        return sort(collect(SS), by = by)
    else
        Λ = [eigvals(S, S.C.dA * [k...], by = nothing) for k in k_iter]
        return Λ = sort(vcat(Λ...), by = by)
    end
end

"""
    eigen(S::X, k; by = abs) where {X<:Union{CrystalOperator,OperatorComposition}}

Computes the eigenvalues and eigenvectors of the symbol of a CrystalOperator/OperatorComposition wrt frequency k.

The eigenvalues are sorted by `by`.
"""
function eigen(S::X, k; by = abs) where {X<:Union{CrystalOperator,OperatorComposition}}
    symb = symbol(S, k)
    ev = LinearAlgebra.eigen(ComplexF64.(symb))
    p = sortperm(ev.values, by = by)
    return LinearAlgebra.Eigen(ev.values[p], ev.vectors[:, p])
end

"""
    eigvals_df(S::X; N = 20, by = abs) where {X<:Union{CrystalOperator,OperatorComposition}}

Returns a dataframe with the eigenvalues of the symbol of a CrystalOperator/OperatorComposition.
The Frequency space is divided into N^dim equidistant (unique) points frequencies k.

"""
function eigvals_df(S::X; N = 20, by = abs) where {X<:Union{CrystalOperator,OperatorComposition}}
    k_iter = Iterators.product([
        range(0, stop = 1, length = N + 1)[1:end-1] for _ = 1:S.C.dim
    ]...)
    kv = reshape([[x...] for x in collect(k_iter)], :)

    dAkv = [S.C.dA * x for x in kv]

    Λ = [eigvals(S, k, by = by) for k in dAkv]

    df = DataFrame(k = kv, dAk = dAkv, Λ = Λ)
end
