using alfa
using Test

# test if show functions throw an error.
@testset "show.jl" begin
    @test sprint(show, MIME("text/plain"), alfa.Crystal()) isa String
    @test sprint(show, MIME("text/plain"), alfa.Multiplier()) isa String
    @test sprint(show, MIME("text/plain"), alfa.gallery.Laplace2D()) isa String
end
