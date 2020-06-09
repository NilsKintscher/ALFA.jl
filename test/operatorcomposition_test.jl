using ALFA
using Test
using LinearAlgebra

@testset "operatorcomposition.jl" begin
        for T in [Float64, Rational{BigInt}]
                 L = ALFA.gallery.Laplace(T=T)
                 LL =L*L
                 f = :($L*$L)
                 oc = ALFA.OperatorComposition(f)
                 k = randn(2)
                 @test ALFA.symbol(oc,k) ≈ ALFA.symbol(LL,k)
                 @test ALFA.eigvals(oc,k) ≈ ALFA.eigvals(LL,k)

                 f = :(I-$L*inv($L) -I +adjoint($L) + I - transpose($L) + I  - I)
                 oc = ALFA.OperatorComposition(f)
                 @test isapprox(norm(ALFA.eigvals(oc,k)),0, atol=1e-14)

        end
end
