using alfa
using Test

@testset "crystal.jl" begin
    @test_throws Exception alfa.Crystal([1 0; 0 1], [1 2 3])
    @test_throws Exception alfa.Crystal([1 0; 0 1], [1 2], [1 2 3])
    @test_throws Exception alfa.Crystal([1 0; 0 1], [1 2 + 1im])

    C = alfa.Crystal()
    @test isa(C, alfa.Crystal) == true
    C = alfa.Crystal([1], [1; 2], [1; 2])
    @test isa(C, alfa.Crystal) == true
    C = alfa.Crystal(nothing, nothing, nothing)
    @test isa(C, alfa.Crystal) == true
    C = alfa.Crystal([1 0; 0 2], [1 2; 3 4; 5 6], [1 2; 3 4; 5.5 6; 7 8])
    @test isa(C, alfa.Crystal) == true

    #getproperty tests
    @test C.size_domain == 3
    @test C.size_codomain == 4
    @test C.dim == 2
    @test C.A == C.L.A
end
