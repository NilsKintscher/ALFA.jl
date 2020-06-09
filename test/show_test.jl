using ALFA
using Test

# test if show functions throw an error.
@testset "show.jl" begin
    @test sprint(show, MIME("text/plain"), ALFA.Crystal{2,Float64}()) isa String
    @test sprint(show, MIME("text/plain"), ALFA.Multiplier()) isa String
    @test sprint(show, MIME("text/plain"), ALFA.gallery.Laplace(N=2,T=Float64)) isa String
end
