# ALFA.jl : Automated Local Fourier Analysis
<!---
tokens which are not needed right now.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://NilsKintscher.github.io/alfa.jl/stable)
[![Coveralls](https://coveralls.io/repos/github/NilsKintscher/alfa.jl/badge.svg?branch=master)](https://coveralls.io/github/NilsKintscher/alfa.jl?branch=master)
-->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://NilsKintscher.github.io/alfa.jl/dev)
[![Build Status](https://travis-ci.com/NilsKintscher/alfa.jl.svg?branch=master)](https://travis-ci.com/NilsKintscher/alfa.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/NilsKintscher/alfa.jl?svg=true)](https://ci.appveyor.com/project/NilsKintscher/alfa-jl)
[![Codecov](https://codecov.io/gh/NilsKintscher/alfa.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/NilsKintscher/alfa.jl)


This is a Julia implementation of the framework described in

[1] Kahl, K., Kintscher, N. Automated local Fourier analysis (aLFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>

The main purpose of this framework is to enable the reliable and easy-to-use analysis of complex methods on repetitive structures, e.g.,  multigrid methods with complex overlapping block smoothers.

Throughout this framework we refer to definitions, theorem, lemma and algorithms of [1].

Important remarks:
- I recommend to read [1] completely, before digging into this framework. The examples in the paper can be found in the documentation of this framework.
- Some unit tests of this framework are done with randomly generated crystal operators. Sometimes for `N` > 2 these tests fail in case the datatype used for the lattice basis and structure elements is  `T::Float64`. Unfortunately this cannot be avoided when using Float64. Thus, the tests don't throw an error if $95\%$ pass. However, For most real applications this module should still work without any internal errors. When in doubt, I recommend using `T::Rational{BigInt}` when possible as the implemented algorithms seem to be reliable. The downside is the slower runtime in comparison to `T::Float64`.
