const alfa_rtol = 1e-4
const alfa_atol = 1e-7

function islessapprox(A::AbstractVector, B::AbstractVector)
    for (a, b) in zip(A, B)
        if !isapprox(a, b, rtol = alfa_rtol, atol = alfa_atol)
            return isless(a, b)
        end
    end
    return isless(length(A), length(B))
end
