using ALFA
using Test
import DataStructures: SortedSet
using StaticArrays
using DataFrames
using LinearAlgebra
const atol = 1e-14


@testset "crystaloperator.jl" begin
    for T in [Float64, Rational{BigInt}]
        @test isa(ALFA.CrystalOperator{2,T}(), ALFA.CrystalOperator) == true

        C = ALFA.Crystal{2,T}(
            [1 0; 0 1],
            [0 0; 0.25 0.25; 0.5 0.5],
            [0 0; 0.25 0.25],
        )
        @test isa(C, ALFA.Crystal) == true

        CO = ALFA.CrystalOperator{2,T}(C)
        @test isa(CO, ALFA.CrystalOperator) == true
        #
        m_wrong_matsize1 = ALFA.Multiplier([1 1], [1 2; 3 4])
        m_wrong_matsize2 = ALFA.Multiplier([1 1], [1 2 3 4; 5 6 7 8])
        m_wrong_dim = ALFA.Multiplier([1 2 3], [1 2 3; 4 5 6])
        @test_throws AssertionError push!(CO, m_wrong_matsize1)
        @test_throws AssertionError push!(CO, m_wrong_matsize2)
        @test_throws AssertionError push!(CO, m_wrong_dim)


        m12 = ALFA.Multiplier([1 2], [1 2 3; 4 5 6])
        @test isa(push!(CO, m12), ALFA.CrystalOperator) == true

        # Define a sorted set which works fine.
        m12 = ALFA.Multiplier([1 2], [1 2 3; 4 5 6])
        m13 = ALFA.Multiplier([1 3], [1 2 3; 4 5 6])


        mymults1 = SortedSet{ALFA.Multiplier}([m13, m12])

        mymults2 = SortedSet{ALFA.Multiplier}()
        push!(mymults2, m12)
        push!(mymults2, m13)

        @test collect(mymults1) == collect(mymults2)
        @test isa(
            ALFA.CrystalOperator{2,T}(C, mymults2),
            ALFA.CrystalOperator,
        ) == true

        mymults3 = SortedSet{ALFA.Multiplier}([m13, m_wrong_dim])
        @test_throws AssertionError ALFA.CrystalOperator{2,T}(C, mymults3)
        mymults3 = SortedSet{ALFA.Multiplier}([m13, m_wrong_matsize1])
        @test_throws AssertionError ALFA.CrystalOperator{2,T}(C, mymults3)
        mymults3 = SortedSet{ALFA.Multiplier}([m13, m_wrong_matsize2])
        @test_throws AssertionError ALFA.CrystalOperator{2,T}(C, mymults3)

        @test isa(ALFA.gallery.Laplace(N=2, T=T), ALFA.CrystalOperator)
        O = ALFA.gallery.Laplace(N=2, T=T)
        @test isa(O, ALFA.CrystalOperator)
        ## do stuff with O.

        # @test show(O.M) == nothing
        # @test show(O) == nothing


        # test find_multiplier
        @test ALFA.find_multiplier(O, [12, 12]) == nothing
        @test ALFA.find_multiplier(O, [0, 0]) isa ALFA.Multiplier


        # test wrtLattice
        S = ALFA.gallery.Laplace(N=2, T=T)


        S2 = ALFA.wrtLattice(S, 2 * S.C.L.A)
        mult = ALFA.Multiplier([0, 0], [-4 1 1 0; 1 -4 0 1; 1 0 -4 1; 0 1 1 -4])
        mult_orig = deepcopy(mult)
        @test ALFA.find_multiplier(push!(S2, mult), [0, 0]).mat == mult_orig.mat
        @test ALFA.find_multiplier(push!(S2, mult, true), [0, 0]).mat ==
              2 * mult_orig.mat
        # @test ALFA.find_multiplier(push!(O, mult), [0, 0]).mat == [-4]
        # @test ALFA.find_multiplier(push!(O, mult, true), [0, 0]).mat == [-8]


        # check of ==
        @test S == ALFA.wrtLattice(S, S.C.L)
        S2 = ALFA.wrtLattice(S, 2 * S.C.L.A)
        @test S2 == deepcopy(S2)
        S3 = deepcopy(S2)
        push!(S3, mult, true)
        @test S2 != S3
        S3 = deepcopy(S2)
        push!(
            S3,
            ALFA.Multiplier(
                [1000, 100],
                [-4 1 1 5; 2 -4 0 1; 1 0 -4 1; 8 1 1 -4],
            ),
            true,
        )
        @test S2 != S3
        @test S != S2

        ##test symbol
        R = ALFA.gallery.fw_restriction(N=2, T=T)
        @test ALFA.symbol(ALFA.CrystalOperator{2,T}(R.C), rand(2)) == [0 0 0 0]

        S2 = ALFA.wrtLattice(S, 2 * S.C.L.A)
        x = [0.5969463404190447, 0.9794807414680586]
        sym = Complex{Float64}[
            -4.0 + 0.0im 1.3452758388167156 - 0.938501249402159im 1.966939802419922 +
                                                                  0.25500474210516544im 0.0 +
                                                                                        0.0im
            1.3452758388167156 + 0.938501249402159im -4.0 + 0.0im 0.0 + 0.0im 1.966939802419922 +
                                                                        0.25500474210516544im
            1.966939802419922 - 0.25500474210516544im 0.0 + 0.0im -4.0 + 0.0im 1.3452758388167156 -
                                                                        0.938501249402159im
            0.0 + 0.0im 1.966939802419922 - 0.25500474210516544im 1.3452758388167156 +
                                                                  0.938501249402159im -4.0 +
                                                                                      0.0im
        ]
        @test ALFA.symbol(S2, x) isa Matrix

        ev = [
            -0.3763088603685248,
            -3.656889100285505,
            -4.343110899714498,
            -7.623691139631474,
        ]
        @test isapprox(ev, ALFA.eigvals(S2, x), atol = atol)
        @test isapprox(
            sort(ev, by = real),
            sort(ALFA.eigvals(S2, x, by = nothing), by = real),
            atol = atol,
        )

        evN2 = [
            1.3322676295501878e-15,
            -2.0,
            -2.0,
            -2.0,
            -2.0,
            -3.9999999999999973,
            -3.999999999999998,
            -3.9999999999999996,
            -4.0,
            -4.0,
            -4.0,
            -6.0,
            -6.0,
            -6.0,
            -6.0,
            -7.999999999999998,
        ]
        @test isapprox(evN2, ALFA.eigvals(S2, N = 2), atol = atol)
        evN2unique =
            Complex[0.0+0.0im, -2.0+0.0im, -4.0+0.0im, -6.0+0.0im, -8.0+0.0im]
        @test isapprox(
            evN2unique,
            ALFA.eigvals(S2, N = 2, unique = true),
            atol = atol,
        )

        ev = ALFA.eigen(S2,x, by=abs)
        sym = ALFA.symbol(S2,x)
        for i = 1:4
            @test sym*ev.vectors[:,i] ≈ ev.values[i]*ev.vectors[:,i]
        end


        k = Array{Float64,1}[[0.0, 0.0], [0.5, 0.0], [0.0, 0.5], [0.5, 0.5]]
        dAk = StaticArrays.SArray{Tuple{2},Float64,1,2}[
            [0.0, 0.0],
            [0.25, 0.0],
            [0.0, 0.25],
            [0.25, 0.25],
        ]
        Λ = Array{Float64,1}[
            [
                1.3322676295501878e-15,
                -3.9999999999999973,
                -4.0,
                -7.999999999999998,
            ],
            [-2.0, -2.0, -6.0, -6.0],
            [-2.0, -2.0, -6.0, -6.0],
            [-3.999999999999998, -3.9999999999999996, -4.0, -4.0],
        ]
        df = DataFrame(k = k, dAk = dAk, Λ = Λ)
        @test isapprox(ALFA.eigvals_df(S2, N = 2).k, k, atol = atol)
        @test isapprox(ALFA.eigvals_df(S2, N = 2).dAk, dAk, atol = atol)
        @test isapprox(ALFA.eigvals_df(S2, N = 2).Λ, Λ, atol = atol)


        for jj = 1:10
            S = rand(ALFA.CrystalOperator{2,T}, single_domain = true)
            Id = ALFA.CrystalOperator(S.C, 2 * I)
            @test ALFA.IsApproxEquivalent((S + Id) - Id, S)
            @test ALFA.IsApproxEquivalent((S + I) - I, S)
            @test !ALFA.IsApproxEquivalent(S + I, S)
            @test !ALFA.IsApproxEquivalent(S + Id, S)
            @test (1 * S) == (S * 1) == S
            @test ALFA.IsApproxEquivalent(2 * S, S * 2)
            @test ALFA.IsApproxEquivalent((S * 4) / 2, 2 * S)
            @test ALFA.IsApproxEquivalent(S^1, S)
            @test S^2 == S * S
        end
    end
end
