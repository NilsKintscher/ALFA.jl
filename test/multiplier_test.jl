using alfa
using Test

@testset "multiplier.jl" begin
    @test isa(alfa.Multiplier(), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1 2]), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1, 2]), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1, 2],  [0]), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1, 2],  [0;1;2]), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1, 2],  [0 1; 1 2]), alfa.Multiplier) == true

    #lexicographically ordered. testing isless and isequal
    m013 = alfa.Multiplier([0 1 4])
    m113 = alfa.Multiplier([1 1 3])
    m123 = alfa.Multiplier([1 2 3])

    @test (m113 < m113) == false
    @test (m113 <= m113) == true
    @test (m113 < m123) == true
    @test (m013 < m113) == true
end
