module alfa

using LinearAlgebra
using StaticArrays
import AbstractAlgebra
import DataStructures: SortedSet, SortedDict
using RecipesBase
using ColorSchemes

const alfa_rtol = 1e-4
const alfa_atol = 1e-7


include("abstractalgebra_wrapper.jl")
include("lattice.jl")
include("crystal.jl")
include("multiplier.jl")
include("crystaloperator.jl")
include("operatorcomposition.jl")
include("gallery.jl")
include("plotrecipes.jl")

end # module
