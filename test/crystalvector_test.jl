using ALFA
using Test
using StaticArrays
using LinearAlgebra
using Random

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

                CV1 = ALFA.CrystalVector(CT, (x,y) -> rand(1:20, x,y))
                CV1C2 = ALFA.wrtLattice(CV1, C2)

                @test ALFA.IsApproxEquivalent(CV1, CV1C2)

                CV1_normCoords = ALFA.ShiftCoordsIntoStandardCell(CV1)
                s_fractional = copy(CV1_normCoords.CT.coords)
                s_cartesian =  [CV1.CT.C.A*x for x in s_fractional]

                CV1_manual_change = ALFA.ChangeTorusCoords(CV1, s_fractional)

                @test ALFA.IsApproxEquivalent(CV1_normCoords, CV1_manual_change)

                CV1_manual_change_cartesian = ALFA.ChangeTorusCoords(CV1, s_cartesian, fractional=false)
                @test ALFA.IsApproxEquivalent(CV1_normCoords, CV1_manual_change_cartesian)

                CV1_manual_change_from_normal = ALFA.ChangeTorusCoords(CV1_normCoords, CV1.CT.coords)

                @test ALFA.IsApproxEquivalent(CV1, CV1_manual_change_from_normal)

                s_fractional_perm = s_fractional[Random.randperm(length(s_fractional))]

                # random
                s_fractional_perm_shift = [x.+1 for x in s_fractional_perm]

                CV1_Change_both = ALFA.ChangeTorusCoords(CV1, s_fractional_perm_shift)

                @test ALFA.IsApproxEquivalent(CV1, CV1_Change_both)


            end
        end
    end
end
