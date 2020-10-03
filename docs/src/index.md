# ALFA.jl : Automated Local Fourier Analysis


[1] Kahl, K., Kintscher, N. Automated local Fourier analysis (aLFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>

The main purpose of this framework is to enable the reliable and easy-to-use analysis of complex methods on repetitive structures, e.g.,  multigrid methods with complex overlapping block smoothers.

Throughout this framework we refer to definitions, theorems, lemmata and algorithms in [1].


## Installation

As this package is registered with the [General Registry](https://github.com/JuliaRegistries/General), it can be installed via Julia's Package Manager as follows:

```
using Pkg
Pkg.add("ALFA")
```

## Some remarks:
- If you only want to use this framework to analyze an operator or a method/composition of operators, you may try to proceed in the same way as in the examples: The examples in [1] can be found in the documentation of this framework.
- If you are interested in the core algorithms and want to understand this framework completely, I recommend to read [1] completely before digging into the source code of this repository.
- Several unit tests of this framework are done with randomly generated crystal operators which can lead to very ill-conditioned problems. Thus, these tests may fail in case the datatype used for the lattice basis and structure elements is  `T::Float64`. As this cannot be avoided when using floating-point arithmetic, the tests don't throw an error if at least $95\%$ pass. From my own experience I can say that this is not a problem in actual applications. Nevertheless, there exist rational versions ( `T::Rational{BigInt}`) of these algorithms which are very reliable, but also slower.
