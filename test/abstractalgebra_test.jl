using alfa
using AbstractAlgebra
using Test


@testset "abstractalgebra_wrapper.jl" begin
    A = [70 1 21; 58 92 28; 64 94 8]
    AZZ = AbstractAlgebra.ZZ[70 1 21; 58 92 28; 64 94 8]
    #convert tests
    @test convert(AbstractAlgebra.Generic.Mat{BigInt}, A) == AZZ
    @test convert(Matrix{Int}, AZZ) == A

    #functions test
    H = [1 0 21406; 0 2 33146; 0 0 70274]
    @test H == alfa.hnf(A)
    S = [1 0 0; 0 2 0; 0 0 70274]
    U = [1 0 0; 244 3647 -3572; 218 3258 -3191]
    V = [0 -1 39332; 1 70 -2753261; 0 0 1]
    @test (S, U, V) == alfa.snf_with_transform(A)
end
