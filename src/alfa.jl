module alfa

using LinearAlgebra
import AbstractAlgebra
import DataStructures: SortedSet, SortedDict
using RecipesBase
using ColorSchemes

include("abstractalgebra_wrapper.jl")
include("lattice.jl")
include("crystal.jl")
include("multiplier.jl")
include("crystaloperator.jl")
include("operatorcomposition.jl")
include("gallery.jl")
include("plotrecipes.jl")

end # module
