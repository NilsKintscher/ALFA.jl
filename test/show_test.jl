using alfa
using Test

@testset "show.jl" begin
    @test sprint(show, MIME("text/plain"), alfa.Crystal()) == "Lattice Basis: alfa.Lattice{2}([1.0 0.0; 0.0 1.0])\nDomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:\n [0.0, 0.0]\nCodomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:\n [0.0, 0.0]"
    @test sprint(show, MIME("text/plain"), alfa.gallery()) ==  "Lattice Basis: alfa.Lattice{2}([1.0 0.0; 0.0 1.0])\nDomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:\n [0.0, 0.0]\nCodomain: 1-element Array{StaticArrays.SArray{Tuple{2},Float64,1,2},1}:\n [0.0, 0.0]\nMultiplier: 5-element Array{alfa.Multiplier,1}:\n alfa.Multiplier{2}([-1, 0], [1])\n alfa.Multiplier{2}([0, -1], [1])\n alfa.Multiplier{2}([0, 0], [-4])\n alfa.Multiplier{2}([0, 1], [1]) \n alfa.Multiplier{2}([1, 0], [1])"
end
