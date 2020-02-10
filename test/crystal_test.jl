using alfa
using Test

@testset "crystal.jl" begin
    @test_throws AssertionError alfa.Crystal([1 0; 0 1], [1 2 3])

    C = alfa.Crystal()
    @test isa(C, alfa.Crystal) == true
    C = alfa.Crystal([1], nothing, [1; 2])
    @test isa(C, alfa.Crystal) == true
    C = alfa.Crystal(nothing, nothing, nothing)
    @test isa(C, alfa.Crystal) == true
    C = alfa.Crystal([1 0; 0 2], [1 2; 3 4; 5 6], nothing)
    @test isa(C, alfa.Crystal) == true
end
