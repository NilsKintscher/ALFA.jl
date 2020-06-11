# Example 3: Half-hybrid smoother for the curl-curl equation

In here we use this framework to analyze the half hybrid smoother for the curl-curl equation as described in [2] <https://epubs.siam.org/doi/abs/10.1137/S0036142997326203>, and analyzed in [3] <https://epubs.siam.org/doi/abs/10.1137/070679119>.

It corresponds to section 6.2 of [1] Kahl, K., Kintscher, N. Automated local Fourier analysis (ALFA). Bit Numer Math (2020). <https://doi.org/10.1007/s10543-019-00797-w>.


## Definition of the operators
```@example curlcurl; continued = true
using ALFA
using LinearAlgebra
using Plots
using SparseArrays
```
#### System & restriction operator

```@example curlcurl
L = ALFA.gallery.curlcurl(.01)
p1 = plot(L, title="L")

R = ALFA.gallery.curlcurl_restriction()
p2 = plot(R, title="R")

P = R'
Lc = R*L*P

(w,h) = p1.attr[:size]
plot(p1, p2, layout=(2,1), size=(w,2h))

```

### Hybrid smoother

#### The Discrete gradient operator

```@example curlcurl
A = [1 0; 0 1]
Domain = [[.5, 0],[0, .5]]
Codomain = [[0,0]]

C = ALFA.Crystal{2,Float64}(A, Domain, Codomain)
Rs = ALFA.CrystalOperator{2,Float64}(C)

push!(Rs, ALFA.Multiplier([0 0], [-1 -1]))
push!(Rs, ALFA.Multiplier([-1 0], [1 0]))
push!(Rs, ALFA.Multiplier([0 -1], [0 1]))

Ps = Rs'
Ls = Rs*L*Ps

p1 = plot(Rs, title="Rs")
p2 = plot(Ps, title="Ps")
p3 = plot(Ls, title="Ls")
plot(p1, p2, p3, layout=(3,1), size=(w,3h))

```

#### Gauss-Seidel on vertices


```@example curlcurl

# We have to change to lexicographic ordering in order to obtain the same operator as described in [3]
S1 = ALFA.CrystalOperatorCopyLowerTriangle(Ls,perm=[2,1])
plot(S1)
```

#### scalar Gauss-Seidel on edges

```@example curlcurl
# change the structure element in order to obtain the same operator as described in [3]
s0_lex = [[.5,0],[1,-.5]]
L_lex = ALFA.ChangeStructureElement(L, s0_lex, s0_lex)
S0 = ALFA.CrystalOperatorCopyLowerTriangle(L_lex, perm=[2,1])

#Now, S0 describes a block-Gauss-Seidel smoother. We need to use the lower triangle of the central multiplier.
m = ALFA.find_multiplier(S0, [0,0])
m.mat = tril(m.mat)

plot(S0)
```



## Spectrum of the smoother

```@example curlcurl
f_edge = :(I-pinv($S0)*$L)
f_node = :(I-$Ps*pinv($S1)*$Rs*$L)

oc_s = ALFA.OperatorComposition(f_edge*f_node)
plotspectrum(oc_s, N=30)
```


## Spectrum of the twogrid method

```@example curlcurl
f_cgc = :(I-$P*inv($Lc)*$R*$L)

oc_tg = ALFA.OperatorComposition(f_node*f_edge*f_cgc*f_node*f_edge)
plotspectrum(oc_tg, N=30)
```


## Double check result with a twogrid implementation

```@example curlcurl
wrtL = ALFA.Lattice{2,Float64}(30*[1.0 0; 0 1.0])
Lm = SparseMatrixCSC{Float64,Int}(ALFA.construct_matrix(L,wrtL))
Rm = SparseMatrixCSC{Float64,Int}(ALFA.construct_matrix(R,wrtL))
Pm = SparseMatrixCSC{Float64,Int}(ALFA.construct_matrix(P,wrtL))
Lcm = SparseMatrixCSC{Float64,Int}(ALFA.construct_matrix(Lc,wrtL))
S0m = SparseMatrixCSC{Float64,Int}(ALFA.construct_matrix(S0,wrtL))
S1m = SparseMatrixCSC{Float64,Int}(ALFA.construct_matrix(S1,wrtL))
Rsm = SparseMatrixCSC{Float64,Int}(ALFA.construct_matrix(Rs,wrtL))
Psm = SparseMatrixCSC{Float64,Int}(ALFA.construct_matrix(Ps,wrtL))

nothing # hide
```


### Definition of the application of the smoother and coarse grid correction

```@example curlcurl
S0mLU = lu(S0m)
S1mLU = lu(S1m)
LcmLU = lu(Lcm)


function smooth_edge(b, x)
    r = b-Lm*x
    x += (S0mLU\r)
    return x
end


function smooth_node(b, x)
    r = b-Lm*x
    rc = Rsm*r
    xc = S1mLU\rc
    x += Psm*xc
    return x
end

function cgc(b, x)
    r = b-Lm*x
    rc = Rm*r
    xc = LcmLU\rc
    x += Pm*xc
    return x
end
nothing # hide
```

### Testrun of  the twogrid method.


```@example curlcurl

# initialize rhs and initial guess
global x, b, casym_vec, resnorm_vec, num_iter
n = size(Lm,1);
x = Lm*rand(n);
x = x/norm(x)
b = 0*Lm*rand(n);

casym_vec = [1.0]
resnorm_vec = [1.0]

num_iter = 0

while resnorm_vec[end] > 1e-200 && num_iter < 175
    global x, b, casym_vec, resnorm_vec, num_iter
    num_iter = num_iter + 1
    x = smooth_edge(b,x)
    x = smooth_node(b,x)
    x = cgc(b,x)
    x = smooth_edge(b,x)
    x = smooth_node(b,x)
    push!(resnorm_vec, norm(b - Lm*x))
    push!(casym_vec, resnorm_vec[end]/resnorm_vec[end-1])
end

# plot convergence behavior.
plot(resnorm_vec, yaxis=:log, xlabel="iteration", ylabel="residual norm", title="measured asymptotic convrate="*string(casym_vec[end])*"\n conv.rate from analysis: "*string(abs(ALFA.eigvals(oc_tg,N=30)[end])))

```
