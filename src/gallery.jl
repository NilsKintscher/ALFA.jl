#
module gallery
using ..alfa
function Laplace2D()
            A = [1 0; 0 1]
            Domain = [0 0]
            Codomain = [0 0]

            C = alfa.Crystal(A, Domain, Codomain)
            L = alfa.CrystalOperator(C)

            push!(L, alfa.Multiplier([0 0], [-4]))
            push!(L, alfa.Multiplier([0 -1], [1]))
            push!(L, alfa.Multiplier([0 1], [1]))
            push!(L, alfa.Multiplier([1 0], [1]))
            push!(L, alfa.Multiplier([-1 0], [1]))
            return L
end

function fw_restriction2D()
            A = 2*[1 0; 0 1]
            Domain = [[0, 0], [0, 1], [1, 0], [1, 1]]
            Codomain = [0, 0]

            C = alfa.Crystal(A, Domain, Codomain)
            L = alfa.CrystalOperator(C)

            push!(L, alfa.Multiplier([0 0], [1 1//2 1//2 1//4]))
            push!(L, alfa.Multiplier([-1 0], [0 0 1//2 1//4]))
            push!(L, alfa.Multiplier([0 -1], [0 1//2 0 1//4]))
            push!(L, alfa.Multiplier([-1 -1], [0 0 0 1//4]))
            return L
end


end
