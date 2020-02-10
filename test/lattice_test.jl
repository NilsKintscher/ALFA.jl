using alfa
using Test

@testset "lattice.jl" begin
    @test_throws DimensionMismatch alfa.Lattice([1 2])
    @test_throws AssertionError alfa.Lattice([1.0 2; 1 2])

    @test isa(alfa.Lattice([5]), alfa.Lattice) == true
    @test isa(alfa.Lattice([5.0]), alfa.Lattice) == true

    Am = [3 1; 5 6]
    A = alfa.Lattice(Am)
    @test isa(A, alfa.Lattice)
    Bm = [7 1; 2 3]
    B = alfa.Lattice(Bm)
    Cm = [-78.0 19.0; -481.0 114.0]
    C = alfa.Lattice(Cm)
    @test lcm(A, B) == C
end
