module alfa

using LinearAlgebra
using StaticArrays
import Nemo
import DataStructures: SortedSet, SortedDict
using RecipesBase
using ColorSchemes
using DataFrames
using MacroTools

# using LazyArrays
#
# import Base.Ordering
# import Base.lt
# import DataStructures.eq
# struct ComplexOrdering <: Ordering
# end
# lt(::ComplexOrdering, a, b) = isless(real(a), real(b)) || (real(a) == real(b) && isless(imag(a), imag(b)))
# eq(::ComplexOrdering, a, b) = a == b


include("utils.jl")
include("nemo_wrapper.jl")
include("lattice.jl")
include("crystal.jl")
include("multiplier.jl")
include("crystaloperator.jl")
include("operatorcomposition.jl")
include("gallery.jl")
include("show.jl")
include("plotrecipes.jl")

end # module
