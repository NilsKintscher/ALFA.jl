
mutable struct OperatorComposition
    C::Crystal
    f::Expr
    OpDict::Dict # Symbol, Any
    sym::Function
end

function OperatorComposition(f::Expr)
    f, dict = wrtSameLatticeAndNormalize(f)
    C = eval(f).C
    f_sym_k = (k; π = π) -> eval(MacroTools.postwalk(f) do x
        if x isa CrystalOperator
            symbol(x, k, π = π)
        else
            x
        end
    end)
    return OperatorComposition(C, f, dict, f_sym_k)
end


function _SimplifyAndApply(f::Expr, g::Function = (x) -> x)
    pulled_reverse_orig = Dict()
    pulled = Dict()
    cnt = 0
    f2 = MacroTools.postwalk(f) do x
        if x isa CrystalOperator
            if x in keys(pulled_reverse_orig)
                d_x = pulled_reverse_orig[x]
            else
                d_x = Symbol(:var_, cnt, "_")
                get!(pulled_reverse_orig, x, d_x)
                gx = g(x)
                get!(pulled, d_x, :($d_x = $gx))
                cnt += 1
            end
            d_x
        else
            x
        end
    end
    out = quote end
    append!(out.args, values(pulled))
    push!(out.args, f2)
    return out, pulled
end

function lcm(f::Expr)
    pulled = Dict()
    cnt = 0
    f2 = MacroTools.postwalk(f) do x
        if x isa CrystalOperator
            d_x = Symbol(:var_, cnt, "_")
            get!(pulled, d_x, x.C.L.A)
            cnt += 1
            nothing
        else
            nothing
        end
    end
    A = lcm(values(pulled)...)
end

function wrtSameLatticeAndNormalize(f::Expr)
    A = lcm(f)
    f2, mydict = _SimplifyAndApply(f, (x) -> wrtLattice(x, A, true))
end
