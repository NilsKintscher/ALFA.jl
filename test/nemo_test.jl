using ALFA
using Nemo
using Test
using StaticArrays
using LinearAlgebra

@testset "nemo_wrapper.jl" begin
    A = [70 1 21; 58 92 28; 64 94 8]
    Astatic = @MMatrix [70 1 21; 58 92 28; 64 94 8]
    AZZ = Nemo.ZZ[70 1 21; 58 92 28; 64 94 8]
    #convert tests
    @test convert(Nemo.fmpz_mat, A) == AZZ
    @test convert(Nemo.fmpz_mat, Astatic) == AZZ
    @test convert(Matrix{BigInt}, AZZ) == A
    @test convert(MMatrix{size(AZZ)...,BigInt}, AZZ) == Astatic
    #functions test
    H = [1 0 0; 0 2 0; 21406 33146 70274]
    Hstatic = @MMatrix [1 0 0; 0 2 0; 21406 33146 70274]

    @test H == ALFA.hnf(A)
    @test Hstatic == ALFA.hnf(Astatic)

    S =  [1 0 0; 0 2 0; 0 0 70274]
    U =  [1 0 0; 244 3647 -3572; 218 3258 -3191]
    V =  [0 -1 39332; 1 70 -2753261; 0 0 1]
    @test (S, U, V) == ALFA.snf_with_transform(A)
    Sstatic = @MMatrix [1 0 0; 0 2 0; 0 0 70274]
    Ustatic = @MMatrix [1 0 0; 244 3647 -3572; 218 3258 -3191]
    Vstatic = @MMatrix [0 -1 39332; 1 70 -2753261; 0 0 1]
    @test (Sstatic, Ustatic, Vstatic) == ALFA.snf_with_transform(Astatic)

    Lstatic = @MMatrix [21 7 -69; 28 -26 34; 8 40 30]
    @test abs(det(Lstatic\Astatic)) ≈ 1
    @test abs(det(Astatic\Lstatic)) ≈ 1

    M = rand(1:100, 10,10)
    Mstatic = convert(MMatrix{size(M)...}, M)

    S = ALFA.lll(M)
    Sstatic = ALFA.lll(Mstatic)
    @test abs(det(S\M)) ≈ 1
    @test abs(det(Sstatic\Mstatic)) ≈ 1
end
