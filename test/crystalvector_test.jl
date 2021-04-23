using ALFA
using Test
using StaticArrays
using LinearAlgebra

@testset "crystalvector.jl" begin
    for T in [Float64, Rational{BigInt}]

        C = ALFA.Crystal{2,T}(
            [1 0; 0 1],
            [0 0; 0.25 0.25; 0.5 0.5],
            [0 0; 0.25 0.25],
        )

        L = ALFA.Lattice{2,T}([2 0; 0 2])
        CT = ALFA.CrystalTorus(C,L)


        for N in 1:3
            for j in 1:10
                C = rand(ALFA.Crystal{N,T}, single_domain=true)
                M1 = rand(-2:3, N, N)
                while(abs(det(M1)) < 1)
                    M1 = rand(-2:3, N, N)
                end
                C2 = C.A*M1

                M2 = rand(-2:3, N, N)
                while(abs(det(M2)) < 1)
                    M2 = rand(-2:3, N, N)
                end
                Z = C2*M2

                CT = ALFA.CrystalTorus(C,Z)

                CV1 = ALFA.CrystalVector(CT, (x) -> rand(1:20, x))
                CV1C2 = ALFA.wrtLattice(CV1, C2)

                @test ALFA.IsApproxEquivalent(CV1, CV1C2)
            end
        end
    end
end
