# alfa.jl : Automated Local Fourier Analysis

This is a julia implementation of the framework described in

[1] Kahl, K., Kintscher, N. Automated local Fourier analysis (aLFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>

The main purpose of this framework is to enable the reliable and easy-to-use analysis of complex methods on repetitive structures, e.g.,  multigrid methods with complex overlapping block smoothers.

Throughout this framework we refer to definitions, theorem, lemma and algorithms of [1].

Important remarks:
- In order to understand this framework, I recommend to read [1] completely. The explained examples in the paper can be found in the documentation of this framework.
- Some unit tests of this framework are done with randomly generated crystal operators. Sometimes for `N` > 2 these tests fail in case the datatype used for the lattice basis and structure elements is  `T::Float64`. Unfortunately this cannot be avoided when using Float64. Thus, the tests don't throw an error if $95\%$ pass. However, For most real applications this module should still work without any internal errors. When in doubt, I recommend using `T::Rational{BigInt}` when possible as the implemented algorithms seem to be reliable. The downside is the slower runtime in comparison to `T::Float64`.
