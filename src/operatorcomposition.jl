
mutable struct OperatorComposition
    f::Expr
    OpDict::Dict # Symbol, Any
    sym::Any
    # list of OperatorCompositions possible..? with reduce()
    function OperatorComposition(f::Expr, g::Dict, h::Any)
        new(f, g, h)
    end
end

function _freevars(literal)
    @show typeof(literal)
    @show literal
    return Symbol[]
end
function _freevars(s::Symbol)
    @show s
    return [s]
end
function _freevars(expr::Expr)
    #esc(expr2)
    if Meta.isexpr(expr, :call)
        return mapfoldl(_freevars, append!, expr.args[2:end], init = Symbol[])
    else
        error("only call expressions allowed!")
    end
end



_freevars2(literal) = literal
function _freevars2(s::Symbol)
    try
        if typeof(eval(s)) <: Function
            return s
        else
            return Meta.parse("DR.x[:$s]")
        end
    catch
        return Meta.parse("DR.x[:$s]")
    end
end
function _freevars2(expr::Expr)
    for (it, s) in enumerate(expr.args)
        expr.args[it] = _freevars2(s)
    end
    return expr
end



macro mymacro4(expr, D)
    expr2 = _freevars2(expr)
    quote
        let DR = Ref($(esc(D)))
            function mysym()
                $expr2
            end
        end
    end
end

function wrtSameLatticeAndNormalize(expr::Expr)
    #@show R
    global vars = _freevars(expr)
    # global x = 2
    # x += 2
    # return x
    # let x = eval(vars[1])
    #     @show x
    # A = nothing
    # for v in vars
    #
    #
    #     if typeof(v) <: CrystalOperator
    #         A = x.C.L.A
    #     else
    #         A = lcm(A, x.C.L.A)
    #     end
    #
    # end
    # for v in vars
    #     v = wrtLattice(v, A)
    # end
    #end
end

macro wrtSameLatticeAndNormalize!(expr)
    vars = (:($(esc(v))) for v in _freevars(expr))
    thearr = []
    quote
        push!(thearr, [$(vars...)])
        #typeof(thedict[1])
        local A = nothing
        for j in thearr
            if typeof(j) <: CrystalOperator
                if A == nothing
                    A = j.C.L.A
                else
                    A = lcm(A, j.C.L.A)
                end
            else
                print(j)
            end
        end
        @show A
        for j in thearr
            if typeof(j) <: CrystalOperator
                j = wrtLattice(j, A)
            else
                j += 5

            end
        end
        return $(QuoteNode(thearr))
        # for v in values($thedict)
        #     $(v) = $(v)*2
        # end
        #thedict
        #return thedict
    end
    #thedict
    #return nothing
end


macro _savefree2(expr)
    fv_pairs =
        (:($(QuoteNode(v)) => $(esc(deepcopy(v)))) for v in _freevars(expr))
    #expr2 = _freevars2(expr)
    thedict = Dict()
    #@show expr2
    @show expr
    f = @mymacro4 expr thedict

    quote
        #thedict2 = Dict($(fv_pairs...))
        push!($(QuoteNode(thedict)), $(fv_pairs...))
        # end
        #     quote

        #@show $(QuoteNode(expr))
        #f = @mymacro4 $(QuoteNode(expr)) thedict
        #
        #
        # f = let DR = Ref($(esc(thedict)))
        #     function mysym()
        #         $expr2
        #     end
        # end

        #f =

        #f = @mymacro4 $(QuoteNode(expr))
        #return :(3+3)
        #OperatorComposition($(QuoteNode(expr)), thedict, $(QuoteNode(f)))
        OperatorComposition(
            $(QuoteNode(expr)),
            $(QuoteNode(thedict)),
            $(QuoteNode(f)),
        )
        #savef(f, thedict)
    end

end



# x = 2
# y = 1
# k = 3.0
#



# function _testsym(s::Symbol)
#     try
#         thetype = typeof(eval(s))
#         @show s
#         @show thetype
#     catch e
#         println("error")
#         @show e
#     end
# end
#
#
# macro _testsymmacro(expr)
#     @show expr
#     quote
#         for (it, s) in enumerate(expr.args)
#             thetype = typeof(eval($(esc(s))))
#             @show thetype
#         end
#     end
#     quote
#         @show typeof($(esc(expr)))
#     end
# end
#
# _freevars3(literal) = literal
# function _freevars3(s::Symbol)
#     try
#         if typeof(eval(s)) <: Function
#             return s
#         else
#             return Meta.parse("DR.x[:$s]")
#         end
#     catch
#         return Meta.parse("DR.x[:$s]")
#     end
# end
#
#
# function _freevars3(expr2::Expr, D::Dict)
#     expr = deepcopy(expr2)
#     for (it, s) in enumerate(expr.args)
#         expr.args[it] = _freevars3(s)
#         # if typeof(s)==Symbol
#         #     expr.args[it] = :(OpDict[$s])
#         #     println(s)
#         #     @show s
#         # end
#     end
#     #myfun = 15
#     myfun = let DR = Ref(D)
#         () -> DR.x[:x] * 2
#     end# eval(expr); end
#     @show expr
#     myfun = let DR = Ref(D)
#         () -> expr
#     end# eval(expr); end
#
#     #myfun
#     return expr
# end
# macro blabla(expr, D)
#     # end)
#     quote
#         let DR = Ref($(esc(D)))
#             #@show DR
#             #() -> eval($(esc(expr)))
#             () -> $(esc(expr))
#         end
#     end
#     #return expr
# end
#
# macro blagen(expr, D)
#     #@show expr
#     expr = _freevars2(expr)
#     @show expr
#
#     quote
#         #expr = 3 * DR.x[:x] + DR.x[:y]
#         # expr2 = _freevars2($(esc(expr)))
#         # @show expr2
#         let DR = Ref($(esc(D))), expr2 = :(3 * DR.x[:x] + DR.x[:y])
#             function mysym()
#                 #DR.x[:x] * 2
#                 $expr2
#             end
#         end
#     end
# end

# macro mymacro1(expr, D)
#     quote
#         let DR = Ref($(esc(D)))
#             function mysym()
#                 3 * DR.x[:x] + DR.x[:y]
#             end
#         end
#     end
# end
#
# macro mymacro2(expr, D)
#     quote
#         let DR = Ref($(esc(D)))
#             function mysym()
#                 $expr
#             end
#         end
#     end
# end
#
# macro mymacro3(expr, D)
#     expr2 = :(3 * DR.x[:x] + DR.x[:y])
#     @show expr2
#     quote
#         let DR = Ref($(esc(D)))
#             function mysym()
#                 $expr2
#             end
#         end
#     end
# end


# #
# # macro mymacro5(expr, D)
# #     quote
# #     expr2 = _freevars2(esc(expr))
# #     quote
# #         let DR = Ref($(esc(D)))
# #             function mysym()
# #                 $expr2
# #             end
# #         end
# #     end
# # end
# # end
# # macro mymacro5(expr, D, pre = "", post = "")
# #     quote
# #         let expr2 = _freevars2($(esc(expr)))#, pre = ", post = post)
# #         let DR = Ref($(esc(D)))
# #             function mysym()
# #                 expr2
# #             end
# #         end
# #     end
# #     end
# # end
#
# macro symbol(expr, D)
#     expr2 = _freevars2(expr)
#     @show expr2
#     quote
#         let DR = Ref($(esc(D)))
#             function mysym(k)
#                 $expr2
#             end
#         end
#     end
# end
#
# function OperatorComposition(f::Expr, OpDict::Dict)
#     sym = @mymacro4 eval(f) OpDict "" ""# "symbol(" ",k)"
#     return OperatorComposition(f, OpDict, sym)
# end
#
# x = 3
# y = 4

#
# struct savef
#     f
#     d
#     function savef(f, d)
#         new(f, d)
#     end
# end
# macro _savefree(expr)
#     fv_pairs =
#         (:($(QuoteNode(v)) => $(esc(deepcopy(v)))) for v in _freevars(expr))
#     #f = @mymacro4
#     oc = :(OperatorComposition($(QuoteNode(expr)), Dict($(fv_pairs...))))
#     #mydict = :(Dict($(fv_pairs...)))
#     #oc = OperatorComposition(expr, Dict{Any,Any}())#, mydict)
#     #oc = OperatorComposition(mydict[1], mydict[2])#, mydict)
#     #println(oc)
#     #return :(OperatorComposition($(QuoteNode(expr), Dict($(fv_pairs...) ))))
#     return oc#OperatorComposition(mydict[1], mydict[2])#eval(QuoteNode(mydict)), eval(QuoteNode(expr))
# end
#
# # macro _saveexpr(expr)
# #     fv_pair = (:(_freevars(expr)))
# #     return fv_pair
# # end
#
#
# macro OperatorComposition(expr) ######### seems to work. Next step: macro for @constructoc oc = expr, suhc that oc.f has oc.D[:x] etc-....
#     fv_pairs =
#         (:($(QuoteNode(v)) => $(esc(deepcopy(v)))) for v in _freevars(expr))
#
#
#     # modify expr
#     expr2 = _freevars2(expr)
#     #
#     mydict = :(OperatorComposition2($(QuoteNode(expr2)), Dict($(fv_pairs...))))
#     #mydict = :(Dict($(fv_pairs...)))
#     #oc = OperatorComposition(expr, Dict{Any,Any}())#, mydict)
#     #oc = OperatorComposition(mydict[1], mydict[2])#, mydict)
#     #println(oc)
#     #return :(OperatorComposition($(QuoteNode(expr), Dict($(fv_pairs...) ))))
#     return mydict#OperatorComposition(mydict[1], mydict[2])#eval(QuoteNode(mydict)), eval(QuoteNode(expr))
# end
#
# function testtest(f)
#     mydict = @_savefree(eval(f))
#     return mydict
# end
# #print(@alfa._savefree(x + 2 * y + k))
# #
# # function mymanipulate(x) = x
# # function mymanipulate(x::Matrix)
# #     return :(Bla[$x])
# # end
# function mymanipulate(x::Symbol)
#     print("\n in symbol \n")
#     print(typeof(eval(x)))
#     return x
# end
#
# ## Constructor of OperatorComposition
# # function OperatorComposition(f::Expr)
# #     #1. Save variables in Dictionary.
# #     OpDict = @_savefree(f)
# #     #2. Manipulate expression
# #
# #     #return OperatorComposition(f, OpDict)
# # end
#
#
#
#
#
# function _ConstructDict!(f::Expr, OpDict::Dict{Symbol,CrystalOperator}) # recursive
#     for (j, key) in enumerate(f.args)
#         # if typeof(key) == Expr
#         #     _ConstructDict!(key, OpDict)
#         # else
#         print("\n eval(key): ")
#         print(Main.eval(key))
#         #     if typeof(key) == Symbol && typeof(eval(key)) == alfa.CrystalOperator
#         #         # check if eval(j) is already present in OpList
#         #         OpName = string(key)
#         #         if haskey(OpDict, key)
#         #             # do not add to dictionary as it already exists.
#         #         else
#         #             push!(OpDict, key => eval(key))
#         #         end
#         #         f.args[j] = :(MyDict[$(QuoteNode(key))])
#         #     end
#         # end
#     end
# end
#
#
#
#
#
# for op in (:*, :+, :-)
#     eval(quote
#         function Base.$op(A::OperatorComposition, B::OperatorComposition)
#             # Consider using deepcopy of the dictionaries here.
#             OpDict = Dict(union(A.OpDict, B.OpDict))
#             f = Meta.parse(string(Expr(:call, $op, A.f, B.f)))
#             return OperatorComposition(f, OpDict)
#         end
#     end)
# end
#
#
#
# macro oc(f::Expr)
#     println(f)
#     return eval(QuoteNode(f)), QuoteNode(f)
# end
#
#
#
# function eigvals(OC::OperatorComposition)
#
# end
#
#
#
# ####
# # given an expression f
# # return dictionary with old variables.
#
# ## macro @alfa.wrtSameLatticeAndNormalize! I-L*2
# #        ; returns L etc noramlized wrtsamelattice.
#
# ## function/macro for symbol: @alfa.symbol I-L [k1,k2]
# ## mapping from oc::OpeatorComposition to symbol(oc, k)?
# ###

function mylcm(expr::Expr)

end

struct ComputationType1 end

# function compute(::ComputationType1, f::Expr)
#     global L
#     #global xx = eval(f.args[1])
#     # for x in f.args
#     #     ex = eval(esc(x))
#     #     if typeof(ex) <: CrystalOperator
#     #         A = lcm(A, ex.C.L.A)
#     #     end
#     # end
#     return L
# end
#
#
#


function getLattice(X::CrystalOperator)
    return X.C.L.A
end
#
# function getLattice(X::Symbol)
#     #return [:(alfa.getLattice($X))]
#     return [X]
# end
#
# function getLattice(X)
#     return []
# end

function getSymbols(X)
    return Symbol[]
end
function getSymbols(X::Symbol)
    return [X]
end
function getSymbols(expr::Expr)
    if Meta.isexpr(expr, :call)
        return mapfoldl(getSymbols, append!, expr.args[2:end], init = Symbol[])
    else
        error("only call expressions allowed!")
    end

end

function getLatticeList(expr::Expr)
    sym_list = unique(getSymbols(expr))
    return [applied(alfa.getLattice, x) for x in sym_list]
    # for (it,x) in enumerate(sym_list)
    #      sym_list[it] = applied(alfa.getLattice, $x) #:(alfa.getLattice($x))
    # end
    #return sym_list
end

# function lcm(a::Vector{Symbol})
#     return [esc(x) for $x in a]#:(lcm($sym_list))
# end
#
# function lcm(expr::Expr)
#     a = getLatticeList(expr)
#     return applied(lcm, a...)
# end
#lcm(expr::Expr) = lcm(eval(expr))
# lcm(expr::Expr, expr::Expr) = lcm(
# function lcm(expr::Expr)
#
# end
#
# function axes(x::alfa.CrystalOperator, j::Int64)
#     return 0
# end
#
# function size(x::alfa.CrystalOperator)
#     return 0
# end
#

# macro ä(ex)
#     LazyArrays.checkex(ex)
#     # Expanding macro here to support, e.g., `@.`
#     esc(:($instantiate2($(LazyArrays.lazy_expr(LazyArrays.macroexpand(__module__, ex))))))
# end
#
# @inline function instantiate2(A::Applied{Style}) where Style
#     #check_applied_axes(A)
#     LazyArrays.Applied{Style}(A.f, map(LazyArrays.instantiate, A.args))
# end
# function LazyArrays.check_applied_axes(A::Applied{<:MatrixFunctionStyle})
#     #length(A.args) == 1 || throw(ArgumentError("MatrixFunctions only defined with 1 arg"))
#     #axes(A.args[1],1) == axes(A.args[1],2) || throw(DimensionMismatch("matrix is not square: dimensions are $axes(A.args[1])"))
# end



struct lattice_lcm_struct
    A
    function lattice_lcm_struct(x)
        new(x)
    end
end

function lattice_lcm(O::CrystalOperator{N,T}) where {N,T}
    lattice_lcm_struct(O.C.L.A)
end

#
lattice_lcm(x::Symbol) = x
#
lattice_lcm(x) = x

# function Base.:*(a::lattice_lcm, b::lattice_lcm)
#     lattice_lcm(Base.lcm(a.A,b.A))
# end
for op in (:*, :+, :-)
    eval(quote
        function Base.$op(a::lattice_lcm_struct, b::lattice_lcm_struct)
            return lattice_lcm_struct(Base.lcm(a.A, b.A))
        end
    end)
end
for op in (:*, :+, :-, :^)
    eval(quote
        function Base.$op(a::lattice_lcm_struct, b)
            return a
        end
    end)
end
for op in (:*, :+, :-)
    eval(quote
        function Base.$op(a, b::lattice_lcm_struct)
            return b
        end
    end)
end

for op in (:transpose, :adjoint, :inv)
    eval(quote
        function Base.$op(a::lattice_lcm_struct)
            return a
        end
    end)
end

wrtLattice(x, y) = x
wrtLattice(x::Expr, y) = x

togglecheck(x) = x
function togglecheck(x::CrystalOperator)
    x._CompatibilityCheckOnly = !(x._CompatibilityCheckOnly)
    return x
end

function lcm(expr::Expr)
    lcmexpr =
        MacroTools.postwalk(x -> x isa Symbol ? x : alfa.lattice_lcm(x), expr)
    A = eval(lcmexpr).A
    return A
end

function mysym(x)
    return
end
function getSymbolExpr(expr::Expr)
    lcmexpr =
        MacroTools.postwalk(x -> x isa Symbol ? x : alfa.lattice_lcm(x), expr)
    A = eval(lcmexpr).A
    rewritexpr =
        MacroTools.postwalk(x -> x isa Symbol ? x : alfa.wrtLattice(x, A), expr)
    checkexpr =
        MacroTools.postwalk(x -> x isa Symbol ? x : alfa.togglecheck(x), rewritexpr)
    try
        eval(checkexpr)
        symexpr = MacroTools.postwalk(
            x -> x isa Symbol ? x : alfa.togglecheck(x),
            checkexpr,
        )
        return symexpr
    catch e
        @show e
    end
end
function symbolfromexpr(expr::Expr, k; π = π)
    symexpr = MacroTools.postwalk(
            x -> x isa CrystalOperator ? alfa.symbol(x, k) : x,
            expr,
        )
        return eval(symexpr)
end

function getsymfun(expr::Expr; π=π)
        return (k;π=π) -> eval(MacroTools.postwalk(
                x -> x isa CrystalOperator ? alfa.symbol(x, k) : x,
                expr,
            ))
end
