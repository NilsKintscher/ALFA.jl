struct Lattice <: AbstractMatrix{Real}
    A::Matrix
    function Lattice(A::Matrix)
        @assert typeof(A) <: Matrix{<:Real} "Matrix entries must be of type Real"
        @assert !isapprox(det(A), 0) "Basis must be nonsingular"
        new(A)
    end
end


function Lattice(mat = nothing)
    if mat == nothing
        mat = Matrix{Float64}(I, 2, 2)
    elseif isa(mat, Vector{T} where {T<:Real}) # if vector convert 1×1 matrix.
        mat = reshape(mat, 1, 1)
    end
    Lattice(mat)
end

Base.size(L::Lattice) = size(L.A)
Base.getindex(L::Lattice, y...) = getindex(L.A, y...)
Base.setindex!(L::Lattice, y...) = setindex!(L.A, y...)

function Base.getproperty(L::Lattice, sym::Symbol)
    if sym == :n || sym == :dim
        size(L.A,1)
    elseif sym == :iA
        inv(L.A)
    elseif sym == :dA
        transpose(inv(L.A))
    else
        # fallback to getfield
        getfield(L, sym)
    end
end


"""
    lcm(A::Matrix,B::Matrix)

Returns the least common multiple of A and B, i.e. a sub-lattice C (C ⊂ A and C ⊂ B) with |det(C)| as small as possible.
"""
function Base.lcm(A::Matrix, B::Matrix, rmax = 50)
    M0 = A \ B
    for r in range(1, stop = rmax)
        rM = r * M0
        rMr = round.(rM)
        if isapprox(rM, rMr)
            # compute SNF
            S, U, V = snf_with_transform(rMr)
            N = diagm(r ./ gcd.(diag(S), r))
            C = B * V * N
            return C
            break
        end
    end
end

function Base.lcm(A::Lattice, B::Lattice, rmax = 50)
    return Lattice(lcm(A.A, B.A))
end
