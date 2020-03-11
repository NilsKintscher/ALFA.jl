using alfa
using Test
import DataStructures: SortedSet
using StaticArrays
using DataFrames
using LinearAlgebra
const atol = 1e-14
@testset "crystaloperator.jl" begin
    @test isa(alfa.CrystalOperator(), alfa.CrystalOperator) == true

    C = alfa.Crystal([1 0; 0 1], [0 0; 0.25 0.25; 0.5 0.5], [0 0; 0.25 0.25])
    @test isa(C, alfa.Crystal) == true

    CO = alfa.CrystalOperator(C)
    @test isa(CO, alfa.CrystalOperator) == true
    #
    m_wrong_matsize1 = alfa.Multiplier([1 1], [1 2; 3 4])
    m_wrong_matsize2 = alfa.Multiplier([1 1], [1 2 3 4; 5 6 7 8])
    m_wrong_dim = alfa.Multiplier([1 2 3], [1 2 3; 4 5 6])
    @test_throws AssertionError push!(CO, m_wrong_matsize1)
    @test_throws AssertionError push!(CO, m_wrong_matsize2)
    @test_throws AssertionError push!(CO, m_wrong_dim)


    m12 = alfa.Multiplier([1 2], [1 2 3; 4 5 6])
    @test isa(push!(CO, m12), alfa.CrystalOperator) == true

    # Define a sorted set which works fine.
    m12 = alfa.Multiplier([1 2], [1 2 3; 4 5 6])
    m13 = alfa.Multiplier([1 3], [1 2 3; 4 5 6])


    mymults1 = SortedSet{alfa.Multiplier}([m13, m12])

    mymults2 = SortedSet{alfa.Multiplier}()
    push!(mymults2, m12)
    push!(mymults2, m13)

    @test collect(mymults1) == collect(mymults2)
    @test isa(alfa.CrystalOperator(C, mymults2), alfa.CrystalOperator) == true

    mymults3 = SortedSet{alfa.Multiplier}([m13, m_wrong_dim])
    @test_throws AssertionError alfa.CrystalOperator(C, mymults3)
    mymults3 = SortedSet{alfa.Multiplier}([m13, m_wrong_matsize1])
    @test_throws AssertionError alfa.CrystalOperator(C, mymults3)
    mymults3 = SortedSet{alfa.Multiplier}([m13, m_wrong_matsize2])
    @test_throws AssertionError alfa.CrystalOperator(C, mymults3)

    @test isa(alfa.gallery.Laplace2D(), alfa.CrystalOperator)
    O = alfa.gallery.Laplace2D()
    @test isa(O, alfa.CrystalOperator)
    ## do stuff with O.

    # @test show(O.M) == nothing
    # @test show(O) == nothing


    # test find_multiplier
    @test alfa.find_multiplier(O, [12, 12]) == nothing
    @test alfa.find_multiplier(O, [0, 0]) isa alfa.Multiplier


    # test wrtLattice
    S = alfa.gallery.Laplace2D()
    S2 = alfa.wrtLattice(S, 2 * S.C.L.A)
    mult = alfa.Multiplier([0, 0], [-4 1 1 0; 1 -4 0 1; 1 0 -4 1; 0 1 1 -4])
    mult_orig = deepcopy(mult)
    @test alfa.find_multiplier(push!(S2, mult), [0, 0]).mat == mult_orig.mat
    @test alfa.find_multiplier(push!(S2, mult, true), [0, 0]).mat ==
          2 * mult_orig.mat
    # @test alfa.find_multiplier(push!(O, mult), [0, 0]).mat == [-4]
    # @test alfa.find_multiplier(push!(O, mult, true), [0, 0]).mat == [-8]


    # check of ==
    @test S == alfa.wrtLattice(S, S.C.L)
    S2 = alfa.wrtLattice(S, 2 * S.C.L.A)
    @test S2 == deepcopy(S2)
    S3 = deepcopy(S2)
    push!(S3, mult, true)
    @test S2 != S3
    S3 = deepcopy(S2)
    push!(
        S3,
        alfa.Multiplier([1000, 100], [-4 1 1 5; 2 -4 0 1; 1 0 -4 1; 8 1 1 -4]),
        true,
    )
    @test S2 != S3
    @test S != S2

    ##test symbol
    R = alfa.gallery.fw_restriction2D()
    @test alfa.symbol(alfa.CrystalOperator(R.C), rand(2)) == [0 0 0 0]

    S2 = alfa.wrtLattice(S, 2 * S.C.L.A)
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
    @test isapprox(sym, alfa.symbol(S2, x), atol = atol)

    ev = [
        -0.3763088603685248,
        -3.656889100285505,
        -4.343110899714498,
        -7.623691139631474,
    ]
    @test isapprox(ev, alfa.eigvals(S2, x), atol = atol)
    @test isapprox(
        sort(ev, by = real),
        sort(alfa.eigvals(S2, x, by = nothing), by = real),
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
    @test isapprox(evN2, alfa.eigvals(S2, N = 2), atol = atol)
    evN2unique =
        Complex[0.0+0.0im, -2.0+0.0im, -4.0+0.0im, -6.0+0.0im, -8.0+0.0im]
    @test isapprox(
        evN2unique,
        alfa.eigvals(S2, N = 2, unique = true),
        atol = atol,
    )

    eigenX = LinearAlgebra.Eigen{
        Complex{Float64},
        Float64,
        Array{Complex{Float64},2},
        Array{Float64,1},
    }(
        [
            -0.3763088603685248,
            -3.656889100285501,
            -4.343110899714497,
            -7.623691139631462,
        ],
        Complex{Float64}[
            -0.44345000006130575 + 0.23098072959800614im -0.4434500000613082 +
                                                         0.23098072959800736im 0.44345000006130586 -
                                                                               0.2309807295980063im 0.4434500000613058 -
                                                                                                    0.23098072959800625im
            -0.49585025491824697 - 0.06428471589351338im 0.49585025491824525 +
                                                         0.06428471589351344im 0.4958502549182475 +
                                                                               0.0642847158935136im -0.4958502549182465 -
                                                                                                    0.06428471589351344im
            -0.4100725299896222 + 0.2860778218385873im -0.4100725299896212 +
                                                       0.2860778218385863im -0.4100725299896234 +
                                                                            0.28607782183858815im -0.4100725299896216 +
                                                                                                  0.2860778218385871im
            -0.5000000000000002 - 0.0im 0.5000000000000007 + 0.0im -0.49999999999999806 -
                                                                   0.0im 0.5000000000000012 -
                                                                         0.0im
        ],
    )

    for i = 1:4
        @test isapprox(
            +eigenX.values[i],
            alfa.eigen(S2, x).values[i],
            atol = atol,
        )
        @test isapprox(
            eigenX.vectors[:, i],
            alfa.eigen(S2, x).vectors[:, i],
            atol = atol,
        ) || isapprox(
            -eigenX.vectors[:, i],
            alfa.eigen(S2, x).vectors[:, i],
            atol = atol,
        )
    end


    k = Array{Float64,1}[[0.0, 0.0], [0.5, 0.0], [0.0, 0.5], [0.5, 0.5]]
    dAk = StaticArrays.SArray{Tuple{2},Float64,1,2}[
        [0.0, 0.0],
        [0.25, 0.0],
        [0.0, 0.25],
        [0.25, 0.25],
    ]
    Λ = Array{Float64,1}[
        [1.3322676295501878e-15, -3.9999999999999973, -4.0, -7.999999999999998],
        [-2.0, -2.0, -6.0, -6.0],
        [-2.0, -2.0, -6.0, -6.0],
        [-3.999999999999998, -3.9999999999999996, -4.0, -4.0],
    ]
    df = DataFrame(k = k, dAk = dAk, Λ = Λ)
    @test isapprox(alfa.compute_spectrum(S2, N = 2).k, k, atol = atol)
    @test isapprox(alfa.compute_spectrum(S2, N = 2).dAk, dAk, atol = atol)
    @test isapprox(alfa.compute_spectrum(S2, N = 2).Λ, Λ, atol = atol)



    S3n = alfa.CrystalOperator(
        alfa.Crystal(
            alfa.Lattice([2.0 0.0; 0.0 2.0]),
            SArray{Tuple{2},Float64,1,2}[
                [0.0, 0.0],
                [0.0, 1.0],
                [1.0, 0.0],
                [1.0, 1.0],
            ],
            SArray{Tuple{2},Float64,1,2}[
                [0.0, 0.0],
                [0.0, 1.0],
                [1.0, 0.0],
                [1.0, 1.0],
            ],
        ),
        SortedSet(
            alfa.Multiplier[
                alfa.Multiplier([-1, 0], [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]),
                alfa.Multiplier([0, -1], [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]),
                alfa.Multiplier(
                    [0, 0],
                    [-4 1 1 0; 1 -4 0 1; 1 0 -4 1; 0 1 1 -4],
                ),
                alfa.Multiplier([0, 1], [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 1 0]),
                alfa.Multiplier([1, 0], [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0]),
                alfa.Multiplier(
                    [1000, 100],
                    [-4 1 1 5; 1 -4 0 1; 2 0 -4 1; 8 1 1 -4],
                ),
            ],
            Base.Order.ForwardOrdering(),
        ),
        false,
    )

    @test S3n == alfa.normalize(S3)

    @test S3n == alfa.normalize(S3n)

    R = alfa.gallery.fw_restriction2D()
    X = alfa.Lattice([0.0 -2.0; 2.0 2.0])

    Rnew = alfa.CrystalOperator(
        alfa.Crystal(
            alfa.Lattice([0.0 -2.0; 2.0 2.0]),
            SArray{Tuple{2},Float64,1,2}[
                [0.0, 0.0],
                [0.0, 1.0],
                [1.0, 0.0],
                [1.0, 1.0],
            ],
            SArray{Tuple{2},Float64,1,2}[[0.0, 0.0]],
        ),
        SortedSet(
            alfa.Multiplier[
                alfa.Multiplier(
                    [-2, 1],
                    Rational{Int64}[0 // 1 0 // 1 0 // 1 1 // 4],
                ),
                alfa.Multiplier(
                    [-1, 0],
                    Rational{Int64}[0 // 1 1 // 2 0 // 1 1 // 4],
                ),
                alfa.Multiplier(
                    [-1, 1],
                    Rational{Int64}[0 // 1 0 // 1 1 // 2 1 // 4],
                ),
                alfa.Multiplier(
                    [0, 0],
                    Rational{Int64}[1 // 1 1 // 2 1 // 2 1 // 4],
                ),
            ],
            Base.Order.ForwardOrdering(),
        ),
        false,
    )
    @test Rnew == alfa.wrtLattice(R, X)

    Rnewn = alfa.CrystalOperator(
        alfa.Crystal(
            alfa.Lattice([0.0 -2.0; 2.0 2.0]),
            StaticArrays.SArray{Tuple{2},Float64,1,2}[
                [-1.0, 1.0],
                [-1.0, 2.0],
                [0.0, 0.0],
                [0.0, 1.0],
            ],
            StaticArrays.SArray{Tuple{2},Float64,1,2}[[0.0, 0.0]],
        ),
        SortedSet(
            alfa.Multiplier[
                alfa.Multiplier(
                    [-1, 0],
                    Rational{Int64}[1 // 4 1 // 2 0 // 1 1 // 2],
                ),
                alfa.Multiplier(
                    [0, -1],
                    Rational{Int64}[1 // 4 1 // 2 0 // 1 0 // 1],
                ),
                alfa.Multiplier(
                    [0, 0],
                    Rational{Int64}[1 // 4 0 // 1 1 // 1 1 // 2],
                ),
                alfa.Multiplier(
                    [1, -1],
                    Rational{Int64}[1 // 4 0 // 1 0 // 1 0 // 1],
                ),
            ],
            Base.Order.ForwardOrdering(),
        ),
        false,
    )

    Snewn = alfa.CrystalOperator(
        alfa.Crystal(
            alfa.Lattice([0.0 -2.0; 2.0 2.0]),
            SArray{Tuple{2},Float64,1,2}[
                [-1.0, 1.0],
                [-1.0, 2.0],
                [0.0, 0.0],
                [0.0, 1.0],
            ],
            SArray{Tuple{2},Float64,1,2}[
                [-1.0, 1.0],
                [-1.0, 2.0],
                [0.0, 0.0],
                [0.0, 1.0],
            ],
        ),
        SortedSet(
            alfa.Multiplier[
                alfa.Multiplier([-1, 0], [0 1 0 0; 0 0 0 0; 0 1 0 1; 0 0 0 0]),
                alfa.Multiplier([-1, 1], [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0]),
                alfa.Multiplier([0, -1], [0 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 0 0]),
                alfa.Multiplier(
                    [0, 0],
                    [-4 1 0 1; 1 -4 0 0; 0 0 -4 1; 1 0 1 -4],
                ),
                alfa.Multiplier([0, 1], [0 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 0]),
                alfa.Multiplier([1, -1], [0 0 0 0; 0 0 0 0; 0 0 0 0; 1 0 0 0]),
                alfa.Multiplier([1, 0], [0 0 0 0; 1 0 1 0; 0 0 0 0; 0 0 1 0]),
            ],
            Base.Order.ForwardOrdering(),
        ),
        false,
    )

    @test Snewn == alfa.wrtSameLatticeAndNormalize(S, R)[1]
    @test Rnewn == alfa.wrtSameLatticeAndNormalize(S, R)[2]
    #
    @test (Snewn,Snewn) == alfa.wrtSameLatticeAndNormalize(Snewn, Snewn)

    #### Test of computation rules.

    SS = alfa.CrystalOperator(alfa.Crystal(alfa.Lattice([1.0 0.0; 0.0 1.0]), SArray{Tuple{2},Float64,1,2}[[0.0, 0.0]], SArray{Tuple{2},Float64,1,2}[[0.0, 0.0]]), SortedSet(alfa.Multiplier[alfa.Multiplier([-2, 0], [1]), alfa.Multiplier([-1, -1], [2]), alfa.Multiplier([-1, 0], [-8]), alfa.Multiplier([-1, 1], [2]), alfa.Multiplier([0, -2], [1]), alfa.Multiplier([0, -1], [-8]), alfa.Multiplier([0, 0], [20]), alfa.Multiplier([0, 1], [-8]), alfa.Multiplier([0, 2], [1]), alfa.Multiplier([1, -1], [2]), alfa.Multiplier([1, 0], [-8]), alfa.Multiplier([1, 1], [2]), alfa.Multiplier([2, 0], [1])],
       Base.Order.ForwardOrdering()), false)

    @test SS == S*S
    @test S == S^1
    @test SS == S^2
    Sc = deepcopy(S)
    Sc._CompatibilityCheckOnly = true

    ScSc = alfa.CrystalOperator(alfa.Crystal(alfa.Lattice([1.0 0.0; 0.0 1.0]), SArray{Tuple{2},Float64,1,2}[[0.0, 0.0]], SArray{Tuple{2},Float64,1,2}[[0.0, 0.0]]), SortedSet(alfa.Multiplier[],
       Base.Order.ForwardOrdering()), true)

    @test Sc*Sc == ScSc

    @test Sc^1 == Sc
    @test Sc^2 == ScSc


end
