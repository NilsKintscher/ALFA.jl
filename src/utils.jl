const alfa_rtol = 1e-4
const alfa_atol = 1e-7

function cmpapprox(A::AbstractVector, B::AbstractVector)
    for (a, b) in zip(A, B)
        if !isapprox(a, b, rtol = alfa_rtol, atol = alfa_atol)
            return isless(a, b) ? -1 : 1
        end
    end
    return cmpapprox(length(A), length(B))
end

islessapprox(A::AbstractVector, B::AbstractVector) = cmpapprox(A, B) < 0
