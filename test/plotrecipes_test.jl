using alfa
using Test
using Plots

# test if plot functions throw an error.
@testset "plotrecipes.jl" begin
    R = alfa.gallery.fw_restriction2D()

    @test plot(R.C.L) isa Plots.Plot{Plots.GRBackend}
    @test plot(R.C) isa Plots.Plot{Plots.GRBackend}
    @test plot(R) isa Plots.Plot{Plots.GRBackend}

    S = alfa.gallery.Laplace2D()
    S2 = alfa.wrtLattice(S, S.C.L.A*2)
    @test surfacespectrum(S2) isa Plots.Plot{Plots.GRBackend}
    @test surfacenorm(S2) isa Plots.Plot{Plots.GRBackend}

end
