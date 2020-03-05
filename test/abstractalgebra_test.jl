using alfa
using AbstractAlgebra
using Test
using StaticArrays

@testset "abstractalgebra_wrapper.jl" begin
    A = [70 1 21; 58 92 28; 64 94 8]
    Astatic = @MMatrix [70 1 21; 58 92 28; 64 94 8]
    AZZ = AbstractAlgebra.ZZ[70 1 21; 58 92 28; 64 94 8]
    #convert tests
    @test convert(AbstractAlgebra.Generic.Mat{BigInt}, A) == AZZ
    @test convert(AbstractAlgebra.Generic.Mat{BigInt}, Astatic) == AZZ
    @test convert(Matrix{Int}, AZZ) == A
    @test convert(MMatrix{size(AZZ)...,Int}, AZZ) == Astatic
    #functions test
    H = [1 0 0; 0 2 0; 21406 33146 70274]
    Hstatic = @MMatrix [1 0 0; 0 2 0; 21406 33146 70274]
    @test H == alfa.hnf(A)
    @test Hstatic == alfa.hnf(A)
    S = [1 0 0; 0 2 0; 0 0 70274]
    Sstatic = @MMatrix [1 0 0; 0 2 0; 0 0 70274]
    U = [1 0 0; 244 3647 -3572; 218 3258 -3191]
    Ustatic = @MMatrix [1 0 0; 244 3647 -3572; 218 3258 -3191]
    V = [0 -1 39332; 1 70 -2753261; 0 0 1]
    Vstatic = @MMatrix [0 -1 39332; 1 70 -2753261; 0 0 1]
    @test (S, U, V) == alfa.snf_with_transform(A)
    @test (Sstatic, Ustatic, Vstatic) == alfa.snf_with_transform(Astatic)
end
