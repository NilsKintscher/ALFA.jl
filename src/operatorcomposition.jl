
struct OperatorComposition
    f::Expr
    OpDict::Dict{Symbol,CrystalOperator}
    # list of OperatorCompositions possible..? with reduce()
    function OperatorComposition(f::Expr, OpDict::Dict{Symbol,CrystalOperator})
        new(f,OpDict)
    end
end


function _ConstructDict!(f::Expr, OpDict::Dict{Symbol,CrystalOperator}) # recursive
    for (j, key) in enumerate(f.args)
        if typeof(key) == Expr
            _ConstructDict!(key, OpDict)
        else
            if typeof(key) == Symbol && typeof(eval(key)) == CrystalOperator
                # check if eval(j) is already present in OpList
                OpName = string(key)
                if haskey(OpDict, key)
                    # do not add to dictionary as it already exists.
                else
                    push!(OpDict, key => eval(key))
                end
                f.args[j] = :(MyDict[$(QuoteNode(key))])
            end
        end
    end
end


function OperatorComposition(f::Expr)
    OpDict = Dict{Symbol,CrystalOperator}()
    _ConstructDict!(f, OpDict)
    return OperatorComposition(f, OpDict)
end


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
