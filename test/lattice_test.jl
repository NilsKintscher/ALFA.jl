using alfa
using Test

@testset "lattice.jl" begin
    @test_throws DimensionMismatch alfa.Lattice([1 2])
    @test_throws AssertionError alfa.Lattice([1 2; 1 2])
    @test_throws AssertionError alfa.Lattice([1 2; 3 4 + 2im])

    @test isa(alfa.Lattice([5]), alfa.Lattice) == true
    @test isa(alfa.Lattice([5.0]), alfa.Lattice) == true

    Am = [3 1; 5 6]
    A = alfa.Lattice(Am)
    @test isa(A, alfa.Lattice)

    @test size(A) == (2, 2)
    @test getindex(A, 1, 2) == 1
    @test A[1, 2] == 1
    @test setindex!(A, 3, 1) == Am
    @test (A[1, 1] = 3) == 3

    invA = [
        0.4615384615384615 -0.07692307692307693
        -0.3846153846153845 0.23076923076923075
    ]
    invAtranspose = [
        0.4615384615384615 -0.3846153846153845
        -0.07692307692307693 0.23076923076923075
    ]
    @test A.iA ≈ invA
    @test A.dA ≈ invAtranspose

    Bm = [7 1; 2 3]
    B = alfa.Lattice(Bm)
    Cm = [-78.0 19.0; -481.0 114.0]
    C = alfa.Lattice(Cm)
    @test lcm(A, B) == C

    A = alfa.Lattice([1 -1; 8 -5])
    B = alfa.Lattice([-1 -4; 7 -5])
    s =    [0 0; -1 -5; -2 -10; -3 -15; -4 -20; -5 -25; -6 -30; -7 -35; -8 -40; -9 -45; -10 -50]
    @test s == alfa.ElementsInQuotientSpace(A,B)
    s2 = [-4.0 -2.0; -4.0 1.0; -3.0 -3.0; -3.0 0.0; -3.0 3.0; -2.0 -1.0; -2.0 2.0; -2.0 5.0; -1.0 1.0; -1.0 4.0; 0.0 0.0]
    y = [-4.0 2.0; -3.0 1.0; -1.0 0.0; -5.0 3.0; -4.0 2.0; -2.0 1.0; -1.0 0.0; -5.0 3.0; -3.0 2.0; -2.0 1.0; 0.0 -0.0]
    p = [9, 6, 3, 11, 8, 5, 2, 10, 7, 4, 1]
    @test (s2,y,p) == alfa.ShiftIntoUnitCell(s,B)

    # more tests needed.

end
