using alfa
using Test
import DataStructures: SortedSet
using StaticArrays
using DataFrames
using LinearAlgebra
const atol = 1e-14


function testit(func::Function, numtests, N, T)
    local cnt = 0
    local errorcase = nothing
    local case = nothing
    local res
    for i = 1:numtests
        (res, case) = eval(func(N, eval(T)))
        if res == false
            #return res, case, cnt
            #cnt = cnt + 1
            errorcase = case
        else
            cnt = cnt + 1
        end
        print(".")
    end
    return true, errorcase, cnt
end



function AB(N, T)
    O = rand(alfa.CrystalOperator{N,T})
    A = rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    B = rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    try
        (An, Bn) = alfa.wrtSameLatticeAndNormalize(A, B)
        if An.C.Domain ≈ Bn.C.Domain &&
           An.C.Codomain ≈ Bn.C.Codomain && An.C.L ≈ Bn.C.L
            return true, nothing
        else
            return false, "Failed with A= $A, B= $B"
        end
    catch
        return false, "Failed with A= $A, B= $B"
    end

end

function commProp(N, T) #(N, T)
    O = rand(alfa.CrystalOperator{N,T})
    A = rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    B = rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)

    try
        if alfa.IsApproxEquivalent(A + B, B + A)
            return true, nothing
        else
            return false, "Failed with A= $A, B= $B"
        end

    catch
        return false, "Failed with A= $A, B= $B"
    end

end

function distProp1(N, T)

    O = rand(alfa.CrystalOperator{N,T})
    B = O#rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    C = rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    A = rand(O, domain_eq_Acodomain = true)


    try
        if alfa.IsApproxEquivalent(A * (B + C), A * B + A * C)
            return true, nothing
        else
            return false, "Failed with A = $A; B= $B; C= $C"
        end
    catch
        return false, "Failed with A = $A; B= $B; C= $C"
    end
end

function distProp2(N, T)
    O = rand(alfa.CrystalOperator{N,T})
    B = O#rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    C = rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    D = rand(O, codomain_eq_Adomain = true)

    try
        if alfa.IsApproxEquivalent((B + C) * D, B * D + C * D)
            return true, nothing
        else
            return false, "Failed with B: $B, C: $C, D: $D"
        end
    catch
        return false, "Failed with B: $B, C: $C, D: $D"
    end
end

function transposeProp(N, T)
    O = rand(alfa.CrystalOperator{N,T})
    A = rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    B = rand(O, codomain_eq_Adomain = true)
    try
        if alfa.IsApproxEquivalent(
            transpose(A * B),
            transpose(B) * transpose(A),
        )
            return true, nothing
        else
            return false, "Failed with A: $A, B: $B"
        end
    catch
        return false, "Failed with A: $A, B: $B"
    end
end

function adjointProp(N, T)
    O = rand(alfa.CrystalOperator{N,T})
    A = rand(O, domain_eq_Adomain = true, codomain_eq_Acodomain = true)
    B = rand(O, codomain_eq_Adomain = true)
    try
        if alfa.IsApproxEquivalent(adjoint(A * B), adjoint(B) * adjoint(A))
            return true, nothing
        else
            return false, "Failed with A: $A, B: $B"
        end
    catch
        return false, "Failed with A: $A, B: $B"
    end
end
# macro check(func::Symbol, x...)
#     quote
#         #@show func
#         local escf = $(esc(func))
#         local xx = $(esc(x))
#         numtests = xx[1]
#         @show xx[2]
#         @show xx[3]
#         local res, case, cnt = testit(escf, xx...)
#         if case === nothing
#             case = "✓"
#         end
#         case = case === nothing ? "✓" : case
#         println("$(nameof(escf)): $cnt / $numtests successful ; $case")
#         @testset "$(nameof(escf)): $cnt / $numtests successful ; $case" begin
#             @test cnt == numtests
#         end
#     end
# end


macro check(func::Symbol, numtests, N, T)
    quote
        local escf = $(esc(func))
        local escnumtests = $(esc(numtests))
        local escN = $(esc(N))
        local escT = $(esc(T))
        local res, case, cnt = testit(escf, escnumtests, escN, escT)
        if case === nothing
            case = "✓"
        end
        case = case === nothing ? "✓" : case
        println("$(nameof(escf)): $cnt / $escnumtests successful ; $case")
        @testset "$(nameof(escf)): $cnt / $escnumtests successful ; $case" begin
            @test cnt == escnumtests
        end
    end
end

numtests = 20
for N in [1,2,3]
    for T in [Float64,Rational{BigInt}]
        @testset "computation properties, crystaloperator (N,T) = ($N,$T)" begin
            @check AB numtests N T
            @check commProp numtests N T
            @check distProp1 numtests N T
            @check distProp2 numtests N T
            @check transposeProp numtests N T
            @check adjointProp numtests N T
        end
    end
end


# macro asdf(x)
#     quote
#         @show esc($x)
#         #local escf = $(esc(func))
#     end
# end
#
# N = 2
# @asdf(2)
