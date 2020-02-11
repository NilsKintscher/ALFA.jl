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
end
