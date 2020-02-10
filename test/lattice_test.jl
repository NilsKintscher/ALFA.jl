using alfa
using Test


@testset "lattice.jl" begin
    @test_throws DimensionMismatch alfa.Lattice([1 2])
    @test_throws AssertionError alfa.Lattice([1 2; 1 2])

    Am = [3 1; 5 6]
    A = alfa.Lattice(Am)
    Bm = [7 1; 2 3]
    B = alfa.Lattice(Bm)
    Cm = [-78.0 19.0; -481.0 114.0]
    C = alfa.Lattice(Cm)
    @test lcm(A,B) == C
end
