const alfa_rtol = 1e-4
const alfa_atol = 1e-7

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
        size(L.A, 1)
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
        if isapprox(rM, rMr, rtol = alfa_rtol, atol = alfa_atol)
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



function ElementsInQuotientSpace(A::Matrix{T}, B::Matrix; return_diag_hnf::Bool=false, return_fractional::Bool = false) where T
    M = A \ B
    Mr = round.(M)
    @assert isapprox(M, Mr, rtol = alfa_rtol, atol = alfa_atol) "A must be a sublattice of B, i.e. A\\B must be integral."
    H = hnf(Mr)
    dH = diag(H)
    if return_diag_hnf
        return dH
    end
    m = prod(dH)
    print(dH)
    J = Iterators.product([0:x-1 for x in dH]...)


    if return_fractional
        s = Matrix{Int}(undef, length(J), length(dH))
        for (i,x) in enumerate(J)
            s[i,:] .= x
        end
    else
        #
        s = Matrix{T}(undef, length(J), length(dH))
        for (i,x) in enumerate(J)
            s[i,:] .= A*[x...]
        end
    end
    return s
end

function ElementsInQuotientSpace(A::Lattice, B::Lattice)
    return ElementsInQuotientSpace(A.A, B.A)
end

function ShiftIntoUnitCell!(s::Matrix, A::Matrix)
    y = transpose(A \ transpose(s))
    map!(
        x -> isapprox(x, round(x), rtol = alfa_rtol, atol = alfa_atol) ?
             round(x) : floor(x),
        y,
        y,
    ) # remove non-integral part ∈(0,1)
    Ay = transpose(A * transpose(y))
    s = s - Ay
    # sort lexicographically
    p = sortperm(
        collect(eachslice(s, dims = 1)),
        alg = Base.Sort.DEFAULT_STABLE,
    )
    s = s[p, :]
    y = y[p, :]
    return s, y, p # We now have s in the primitive cell of A (A*[0,1)^dim) and s is lexicographically ordered
end

function ShiftIntoUnitCell(s::Matrix, A::Matrix)
    s2 = deepcopy(s)
    return ShiftIntoUnitCell!(s2, A)
end

function ShiftIntoUnitCell(s::Matrix, A::Lattice)
    return ShiftIntoUnitCell(s, A.A)
end

function ShiftIntoUnitCell!(s::Matrix, A::Lattice)
    return ShiftIntoUnitCell!(s, A.A)
end
