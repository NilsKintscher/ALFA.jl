using alfa
using Test

@testset "multiplier.jl" begin
    @test isa(alfa.Multiplier(), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1 2]), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1, 2]), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1, 2], [0]), alfa.Multiplier) == true
    @test isa(alfa.Multiplier([1, 2], [0; 1.1; 2]), alfa.Multiplier) == true
    @test isa(
        alfa.Multiplier([1, 2], [2 + 1im 1; 1 2]),
        alfa.Multiplier,
    ) == true

    m013 = alfa.Multiplier([0 1 4])
    m113 = alfa.Multiplier([1 1 3], [1 2; 3 4; 5 6])

    @test m113 == deepcopy(m113)
    @test m013 != m113

    m123 = alfa.Multiplier([1 2 3])

    # test property functions
    @test m013.dim == m013.n == 3
    @test size(m113) == (3, 2)
    @test size(m113, 1) == 3
    @test size(m113, 2) == 2

    @test m113.size_domain == 2
    @test m113.size_codomain == 3

    #lexicographically ordered. testing isless and isequal
    @test isless(m113, m113) == false
    @test isequal(m113, m113) == true
    @test (m113 < m113) == false
    @test (m113 == m113) == true
    @test (m113 <= m113) == true
    @test (m113 < m123) == true
    @test (m013 < m113 <= m113 < m123) == true
end
