using ALFA
using Test
using Plots

# test if plot functions throw an error.
@testset "plotrecipes.jl" begin
        R = ALFA.gallery.fw_restriction(N = 1, T = Float64)

        @test plot(R.C.L) isa Plots.Plot{Plots.GRBackend}
        @test plot(R.C) isa Plots.Plot{Plots.GRBackend}
        @test plot(R) isa Plots.Plot{Plots.GRBackend}

        S = ALFA.gallery.Laplace(N = 1, T = Float64)
        S2 = ALFA.wrtLattice(S, S.C.L.A * 2)
        @test plotspectrum(S2) isa Plots.Plot{Plots.GRBackend}
        @test plotspectrum(S2, 4) isa Plots.Plot{Plots.GRBackend}





        R = ALFA.gallery.fw_restriction(N = 2, T = Float64)

        @test plot(R.C.L) isa Plots.Plot{Plots.GRBackend}
        @test plot(R.C) isa Plots.Plot{Plots.GRBackend}
        @test plot(R) isa Plots.Plot{Plots.GRBackend}

        S = ALFA.gallery.Laplace(N = 2, T = Float64)
        S2 = ALFA.wrtLattice(S, S.C.L.A * 2)
        @test plotspectrum(S2) isa Plots.Plot{Plots.GRBackend}
        @test plotspectrum(S2, 4) isa Plots.Plot{Plots.GRBackend}
end
