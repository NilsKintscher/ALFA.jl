using SafeTestsets

@safetestset "AbstractAlgebra Wrapper" begin
    include("abstractalgebra_test.jl")
end

@safetestset "Lattice" begin
    include("lattice_test.jl")
end
