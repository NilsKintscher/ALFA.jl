function gallery(s = nothing)
    if s == nothing
        s = "Laplacian2D"
    end
    if s == "Laplacian2D"
        A = [1 0; 0 1]
        Domain = [0 0]
        Codomain = [0 0]

        C = alfa.Crystal(A, Domain, Codomain)
        L = alfa.CrystalOperator(C)

        push!(L, Multiplier([0 0], [-4]))
        push!(L, Multiplier([0 -1], [1]))
        push!(L, Multiplier([0 1], [1]))
        push!(L, Multiplier([1 -1], [1]))
        push!(L, Multiplier([1 1], [1]))
        return L
    end
end
