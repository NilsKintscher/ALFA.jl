using ALFA
using Test
using StaticArrays

@testset "crystaltorus.jl" begin
    for T in [Float64, Rational{BigInt}]

        C = ALFA.Crystal{2,T}(
            [1 0; 0 1],
            [0 0; 0.25 0.25; 0.5 0.5],
            [0 0; 0.25 0.25],
        )

        L = ALFA.Lattice{2,T}([2 0; 0 2])
        CT = ALFA.CrystalTorus(C,L)
        @test isa(CT, ALFA.CrystalTorus) == true

    end
end
