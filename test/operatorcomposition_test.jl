using alfa
using Test
using LinearAlgebra

@testset "operatorcomposition.jl" begin
        for T in [Float64, Rational{BigInt}]
                 L = alfa.gallery.Laplace(T=T)
                 LL =L*L
                 f = :($L*$L)
                 oc = alfa.OperatorComposition(f)
                 k = randn(2)
                 @test alfa.symbol(oc,k) ≈ alfa.symbol(LL,k)
                 @test alfa.eigvals(oc,k) ≈ alfa.eigvals(LL,k)

                 f = :(I-$L*inv($L) -I +adjoint($L) + I - transpose($L) + I  - I)
                 oc = alfa.OperatorComposition(f)
                 @test isapprox(norm(alfa.eigvals(oc,k)),0, atol=1e-20)

        end
end
