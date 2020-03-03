using alfa
using Test

@testset "operatorcomposition.jl" begin

    O1 = alfa.CrystalOperator()

    f = :(1*O1)

    C = alfa.OperatorComposition(f)
    @test isa(C, alfa.OperatorComposition) == true





end
