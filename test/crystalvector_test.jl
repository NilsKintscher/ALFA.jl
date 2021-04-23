using ALFA
using Test
using StaticArrays

@testset "crystalvector.jl" begin
    for T in [Float64, Rational{BigInt}]

        C = ALFA.Crystal{2,T}(
            [1 0; 0 1],
            [0 0; 0.25 0.25; 0.5 0.5],
            [0 0; 0.25 0.25],
        )

        L = ALFA.Lattice{2,T}([2 0; 0 2])
        CT = ALFA.CrystalTorus{2,T}(C,L)


        CV = ALFA.CrystalVector(CT, x -> rand(1:20,x))
        @test isa(CV, ALFA.CrystalVector) == true
    end
end
