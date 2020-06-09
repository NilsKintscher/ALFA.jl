# Example 2: Colored overlapping smoother for graphene

In here we use the ALFA framework to analyze a 4 color overlap smoother for the tight-binding Hamiltonian of graphene.

It corresponds to section 6.1 of [1] Kahl, K., Kintscher, N. Automated local Fourier analysis (ALFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>.

## Definition of the operators
```@example ex2; continued = true
using ALFA
using LinearAlgebra
using Plots
```
#### System operator


```@example ex2
L = ALFA.gallery.graphene_tight_binding()
plot(L)
```

#### Restriction operator of the two-grid method
```@example ex2
R = ALFA.gallery.graphene_dirac_restriction(wl=.25, wlh=.25)
plot(R)
```

#### Prolongation and coarse grid operator

```@example ex2
P = R'
Lc = R*L*P

p1 = plot(P, title="P",  aspect_ratio=:equal)
p2 = plot(Lc, title="Lc",  aspect_ratio=:equal)
(w,h) = p1.attr[:size]
plot(p1, p2, layout=(2,1), size=(w,2h))
```

#### Smoother definition

We first rewrite the system operator $L$ with respect to the translational invariance 2A, where $A=$``L.C.L.A``

```@example ex2
slat = ALFA.Lattice{2,Float64}(2*L.C.L.A)
Ls = ALFA.wrtLattice(L,slat)
plot(Ls)
```

We now construct the $4$ operators used in the four color smoother.

```@example ex2
idx=[2,3,4,5,6,7] #
shifts = [[i,j] for i in [0,1] for j in [0,1]]
S = []
p = []
for (it,x) in enumerate(shifts)
   se = [y + L.C.L.A*x for y in Ls.C.Domain]
   Ls_tmp = ALFA.ChangeStructureElement(Ls, se, se)
   push!(S, ALFA.CrystalOperatorCopyWithMultipliers(Ls_tmp, idx=idx))
   push!(p, plot(S[it], title="S[$it]"))
end
plot(p..., layout=(2,2), size=(1.5w,1.5h))
```

... and normalize the operators S for illustrative purposes:

```@example ex2
p = []
for (j,s) in enumerate(S)
   S[j] = ALFA.normalize(s)
   push!(p, plot(S[j], title="S[$j]"))
end
plot(p..., layout=(2,2), size=(1.5w,1.5h))
```

## Spectrum of the error propagator of the Smoother
Now we construct the error propagator of the smoother and plot its spectrum

```@example ex2
f_s = prod(:(I-0.5*pinv($s)*$L) for s in S)
oc_s = ALFA.OperatorComposition(f_s)
surfacespectrum(oc_s, N=41)
```

## Spectrum of the error propagator of the two-grid method
Finally, we analyze the two-grid error propagator

```@example ex2
f_cgc = :(I-$P*inv($Lc)*$R*$L)
f_tg = f_s*f_cgc*f_s
oc_tg = ALFA.OperatorComposition(f_tg)
surfacespectrum(oc_tg, N=41)
```
