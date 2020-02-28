
struct OperatorComposition
    f::Expr
    OpDict::Dict{Symbol,Any}
    # list of OperatorCompositions possible..? with reduce()
    function OperatorComposition(f::Expr, OpDict::Dict{Symbol,Any})
        new(f, OpDict)
    end
end

#
# x = [2 2]
# y = 1
# k = 3.0
# _freevars(literal) = Symbol[]
# _freevars(s::Symbol) = [s]
# function _freevars(expr::Expr)
#   if Meta.isexpr(expr, :call)
#     return mapfoldl(_freevars, append!, expr.args[2:end], init = Symbol[])
#   else
#     error("only call expressions allowed!")
#   end
# end
#
# macro _savefree(expr)
#   fv_pairs = (:($(QuoteNode(v)) => $(esc(deepcopy(v)))) for v in _freevars(expr))
#   return :(Dict($(fv_pairs...)))
# end
#
# print(@alfa._savefree(x + 2f(y) + k))
#
# function mymanipulate(x) = x
# function mymanipulate(x::Matrix)
#     return :(Bla[$x])
# end
# function mymanipulate(x::Symbol)
#     print("\n in symbol \n")
#     print(typeof(eval(x)))
#     return x
# end
#
# ## Constructor of OperatorComposition
# function OperatorComposition(f::Expr)
#     #1. Save variables in Dictionary.
#     OpDict = @_savefree(f)
#     #2. Manipulate expression
#
#     return OperatorComposition(f, OpDict)
# end
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
#             print("\n eval(key): ")
#             print(Main.eval(key))
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



for op in (:*, :+, :-)
    eval(quote
        function Base.$op(A::OperatorComposition, B::OperatorComposition)
            # Consider using deepcopy of the dictionaries here.
            OpDict = Dict(union(A.OpDict, B.OpDict))
            f = Meta.parse(string(Expr(:call, $op, A.f, B.f)))
            return OperatorComposition(f, OpDict)
        end
    end)
end
