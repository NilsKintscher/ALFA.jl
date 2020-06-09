using ALFA
using Test
import DataStructures: SortedSet
using StaticArrays
using DataFrames
using LinearAlgebra
const atol = 1e-14


@testset "gallery.jl" begin
    T = Float64
    for N in 1:2
    L = ALFA.gallery.Laplace(N=N, T=T)
    @test isa(L, ALFA.CrystalOperator) == true
    R = ALFA.gallery.fw_restriction(N=N, T=T, m=1)
    @test isa(R, ALFA.CrystalOperator) == true
    @test isa(R*L*R', ALFA.CrystalOperator) == true
end
    L = ALFA.gallery.graphene_tight_binding(T=T)
    @test isa(L, ALFA.CrystalOperator) == true
    R = ALFA.gallery.graphene_dirac_restriction(T=T, m=1)
    @test isa(R, ALFA.CrystalOperator) == true
    @test isa(R*L*R', ALFA.CrystalOperator) == true
end
