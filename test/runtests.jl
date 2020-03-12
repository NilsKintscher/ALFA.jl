using SafeTestsets

@safetestset "Nemo Wrapper" begin
    include("nemo_test.jl")
end

@safetestset "Lattice" begin
    include("lattice_test.jl")
end

@safetestset "Crystal" begin
    include("crystal_test.jl")
end

@safetestset "Multiplier" begin
    include("multiplier_test.jl")
end

@safetestset "CrystalOperator" begin
    include("crystaloperator_test.jl")
end

@safetestset "OperatorComposition" begin
    include("operatorcomposition_test.jl")
end


@safetestset "Show" begin
    include("show_test.jl")
end

@safetestset "PlotRecipes" begin
    include("plotrecipes_test.jl")
end
