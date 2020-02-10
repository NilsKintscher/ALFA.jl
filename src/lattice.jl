struct Lattice <: AbstractMatrix{Float64}
    A::Matrix{Float64}
    function Lattice(A::Matrix{Float64})
        @assert !isapprox(det(A), 0)
        new(A)
    end
end

function Lattice(mat = nothing)
    if mat == nothing
        new(Matrix{Float64}(I, 2, 2))
    else
        mat = convert(Matrix{Float64}, mat)
        Lattice(mat)
    end
end

Base.size(L::Lattice) = size(L.A)
Base.getindex(L::Lattice, y...) = getindex(L.A, y...)
Base.setindex!(L::Lattice, y...) = setindex!(L.A, y...)

function Base.getproperty(L::Lattice, sym::Symbol)
    if sym == :n || sym == :dim
        size(L.A)[1]
    elseif sym == :iA
        inv(L.A)
    elseif sym == :dA
        transpose(inv(L.A))
    else
        # fallback to getfield
        getfield(L, sym)
    end
end


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
