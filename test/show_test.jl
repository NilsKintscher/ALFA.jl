using alfa
using Test

# test if show functions throw an error.
@testset "show.jl" begin
    @test sprint(show, MIME("text/plain"), alfa.Crystal{2,Float64}()) isa String
    @test sprint(show, MIME("text/plain"), alfa.Multiplier()) isa String
    @test sprint(show, MIME("text/plain"), alfa.gallery.Laplace(2,Float64)) isa String
end
