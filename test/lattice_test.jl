using alfa
using Test
using StaticArrays
using LinearAlgebra

@testset "lattice.jl" begin

    ## types
    for T in [Float64, Rational{Int}, Rational{BigInt}]
        @test_throws AssertionError alfa.Lattice{1,T}([1 2])
        @test_throws AssertionError alfa.Lattice{2,T}([1 2; 1 2])

        @test isa(alfa.Lattice{1,T}([5]), alfa.Lattice) == true

        #check if the output of print() can be used as a constructor.


        Amstatic = @MMatrix [3 1; 5 6]
        A = alfa.Lattice{2,T}(Amstatic)
        @test isa(A, alfa.Lattice)
        # testing == operator
        @test A == deepcopy(A)
        @test A != alfa.Lattice{1,T}([5])

        @test size(A) == (2, 2)
        @test getindex(A, 1, 2) == 1
        # @test A[1, 2] == 1
        # @test (A[1, 1] = 3) == 3

        @test norm(A.iA * A.A - I) == 0
        @test norm(A.dA * A.A' - I) == 0

        Bm = [7 1; 2 3]
        B = alfa.Lattice{2,T}(Bm)

        @test B.iA * lcm(A, B).A ≈ round.(B.iA * lcm(A, B).A)
        @test A.iA * lcm(A, B).A ≈ round.(A.iA * lcm(A, B).A)

        A = alfa.Lattice{2,T}([1 -1; 8 -5])
        B = alfa.Lattice{2,T}([-1 -4; 7 -5])
        s = [
            SVector{2,T}(x)
            for
            x in eachrow([
                0 0
                -1 -5
                -2 -10
                -3 -15
                -4 -20
                -5 -25
                -6 -30
                -7 -35
                -8 -40
                -9 -45
                -10 -50
            ])
        ]
        sfrac = StaticArrays.SArray{Tuple{2},T,1,2}[
            [0, 0],
            [0, 1],
            [0, 2],
            [0, 3],
            [0, 4],
            [0, 5],
            [0, 6],
            [0, 7],
            [0, 8],
            [0, 9],
            [0, 10],
        ]
        @test s == alfa.ElementsInQuotientSpace(A.A, B.A)
        @test sfrac ==
              alfa.ElementsInQuotientSpace(A.A, B.A, return_fractional = true)
        @test s ==
              alfa.ElementsInQuotientSpace(A.A, B.A, return_diag_hnf = true)[1]
        @test [1, 11] ==
              alfa.ElementsInQuotientSpace(A.A, B.A, return_diag_hnf = true)[2]

        s2 = [
            SVector{2,Float64}(x)
            for
            x in eachrow([
                -4.0 -2.0
                -4.0 1.0
                -3.0 -3.0
                -3.0 0.0
                -3.0 3.0
                -2.0 -1.0
                -2.0 2.0
                -2.0 5.0
                -1.0 1.0
                -1.0 4.0
                0.0 0.0
            ])
        ]
        y = [
            SVector{2,Int}(x)
            for
            x in eachrow([
                -4.0 2.0
                -3.0 1.0
                -1.0 0.0
                -5.0 3.0
                -4.0 2.0
                -2.0 1.0
                -1.0 0.0
                -5.0 3.0
                -3.0 2.0
                -2.0 1.0
                0.0 -0.0
            ])
        ]
        p = [9, 6, 3, 11, 8, 5, 2, 10, 7, 4, 1]
        (s2, y, p) == alfa.ShiftIntoUnitCell(s, B)

        @test norm([floor.(B.iA * x) for x in s2]) == 0
        @test norm([ss + B.A * ys - s[pp] for (ss, ys, pp) in zip(s2, y, p)]) == 0
        # more tests needed.
    end
end
