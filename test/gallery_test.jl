using alfa
using Test
import DataStructures: SortedSet
using StaticArrays
using DataFrames
using LinearAlgebra
const atol = 1e-14


@testset "gallery.jl" begin
    T = Float64
    for N in 1:2
    L = alfa.gallery.Laplace(N=N, T=T)
    @test isa(L, alfa.CrystalOperator) == true
    R = alfa.gallery.fw_restriction(N=N, T=T, m=1)
    @test isa(R, alfa.CrystalOperator) == true
    @test isa(R*L*R', alfa.CrystalOperator) == true
end
    L = alfa.gallery.graphene_tight_binding(T=T)
    @test isa(L, alfa.CrystalOperator) == true
    R = alfa.gallery.graphene_dirac_restriction(T=T, m=1)
    @test isa(R, alfa.CrystalOperator) == true
    @test isa(R*L*R', alfa.CrystalOperator) == true
end
