module alfa

using LinearAlgebra
import AbstractAlgebra
import DataStructures: SortedSet

include("abstractalgebra_wrapper.jl")
include("lattice.jl")
include("crystal.jl")
include("multiplier.jl")
include("crystaloperator.jl")
include("operatorcomposition.jl")

end # module
