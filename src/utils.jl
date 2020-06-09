const ALFA_rtol = 1e-4
const ALFA_atol = 1e-7

function islessapprox(A::AbstractVector, B::AbstractVector)
    for (a, b) in zip(A, B)
        if !isapprox(a, b, rtol = ALFA_rtol, atol = ALFA_atol)
            return isless(a, b)
        end
    end
    return isless(length(A), length(B))
end
