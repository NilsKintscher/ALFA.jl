module alfa

using LinearAlgebra
import AbstractAlgebra

include("abstractalgebra_wrapper.jl")
include("lattice.jl")

"""
    func(x)

Returns double the number `x` plus `1`.
"""
func(x) = 2x + 1

end # module
