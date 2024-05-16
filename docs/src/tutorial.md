# Tutorial

## Installation

At the moment the package is not registered,
so you can install it by running the following command in the Julia REPL:

```julia-repl
pkg> no-you-can-not-install-it-yet
```

## Basic functionalities

To start using QuiverTools in the REPL, one first must import it.

```julia-repl
julia> using QuiverTools
```

Quivers can be built by passing the adjacency matrix to the `Quiver()` constructor:

```julia-repl
julia> Quiver([0 3; 0 0])
Quiver with adjacency matrix [0 3; 0 0]
```

The constructor accepts an optional string for naming the quiver:

```julia-repl
julia> MyQ = Quiver([0 3; 0 0], "My personal quiver")
My personal quiver, with adjacency matrix [0 3; 0 0]
```

QuiverTools has several constructors in place for many common examples:

```julia-repl
julia> mKronecker_quiver(4)
4-Kronecker quiver, with adjacency matrix [0 4; 0 0]

julia> loop_quiver(5)
5-loop quiver, with adjacency matrix [5;;]

julia> subspace_quiver(3)
3-subspace quiver, with adjacency matrix [0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 0 0]

julia> three_vertex_quiver(1, 6, 7)
An acyclic 3-vertex quiver, with adjacency matrix [0 1 6; 0 0 7; 0 0 0]
```

Dimension vectors and stability parameters are represented by `Vector{Int}` objects:

```julia-repl
julia> Q = mKronecker_quiver(3); d = [2,3];

julia> θ = canonical_stability(Q, d)
2-element Vector{Int64}:
  9
 -6

julia> is_coprime(d, θ)
true
```

Here, `is_coprime()` checks if ``d`` is θ-coprime, i.e., if none of the
proper subdimension vectors ``0 \neq d' \nleq d`` satisfies ``\theta \cdot d' = 0``.

The bilinear Euler form relative to a quiver Q of any two vectors
in ``\mathbb{Z}^{Q_0}`` can be computed:

```julia-repl
julia> Q = mKronecker_quiver(3); d = [2,2]; e = [3,4];

julia> Euler_form(Q, d, e)
-10

julia> Euler_form(Q, e, d)
-4
```

This allows to verify whether any given dimension vector belongs to the fundamental domain of the problem.

The fundamental domain is the cone of dimension vectors in ``\mathbb{Z}^{Q_0}``
such that the symmetric Tits form is negative on all the simple roots, i.e.,
for all vertices i,

```math
(s_i, d) := \langle d, s_i\rangle + \langle s_i, d\rangle  \leq 0,
```

where ``s_i`` is the dimension vector with all entries set to ``0`` and the i-th
set to ``1``.

```julia-repl
julia> QuiverTools.in_fundamental_domain(Q, d)
true

julia> QuiverTools.in_fundamental_domain(Q, [1,3])
false
```

## Stable and semistable dimension vectors

One can check if semistable, respectively stable representations
exist for a given dimension vector and stability parameter:

```julia-repl
julia> Q = mKronecker_quiver(3); d = [2,3]; θ = [3,-2];

julia> has_semistables(Q, d, θ)
true

julia> has_stables(Q, d, θ)
true

julia> K2 = mKronecker_quiver(2);

julia> has_stables(K2, [2,2], [1,-1])
false

julia> has_semistables(K2, [2,2], [1,-1])
true
```

One can also determine whether stable representations exist at all
for a given dimension vector by checking if it is a Schur root:

```julia-repl
julia> Q = mKronecker_quiver(3); d = [2, 2];

julia> QuiverTools.is_Schur_root(Q, d)
true

julia> K2 = mKronecker_quiver(2);

julia> QuiverTools.is_Schur_root(K2, d)
false
```

To investigate the Harder-Narasimhan stratification of the parameter space
``\mathrm{R}(Q,\mathbf{d})``, the module provides an implementation based on a recursive algorithm.

```julia-repl
julia> Q = mKronecker_quiver(3); d = [2, 3]; θ = [3, -2];

julia> allHNtypes(Q, d, θ)
8-element Vector{Vector{Vector{Int64}}}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> is_amply_stable(Q, d, θ)
true
```

The method `is_amply_stable()` determines whether the codimension of the θ-semistable locus,
``\mathrm{R}^{\theta-sst}(Q,\mathbf{d})\subset\mathrm{R}(Q,\mathbf{d})``, is at least 2.

The method `all_HN_types()` provides a list of all the Harder-Narasimhan types that appear in the problem.

The method `all_Teleman_bounds()` computes the bounds to apply Teleman quantization on the non-dense strata.
The output is a dictionary whose keys are the HN types and whose values are the weights themselves.

```julia-repl
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> all_Teleman_bounds(Q, d, theta)
Dict{Vector{Vector{Int64}}, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 20
  [[2, 1], [0, 2]]         => 100
  [[1, 0], [1, 2], [0, 1]] => 100
  [[1, 0], [1, 3]]         => 120
  [[1, 0], [1, 1], [0, 2]] => 90
  [[1, 1], [1, 2]]         => 15
  [[2, 0], [0, 3]]         => 90
```

## Verify Teleman inequalities

In the following example, for each ``i,j`` and on each Harder-Narasimhan stratum,
we compute the weight of ``\mathcal{U}_i^\vee \otimes \mathcal{U}_j`` relative to the
1-PS corresponding to the HN stratum. These are then compared to the Teleman bounds.

```julia-repl
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> hn = all_Teleman_bounds(Q,d,theta)
Dict{Vector{Vector{Int64}}, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 20
  [[2, 1], [0, 2]]         => 100
  [[1, 0], [1, 2], [0, 1]] => 100
  [[1, 0], [1, 3]]         => 120
  [[1, 0], [1, 1], [0, 2]] => 90
  [[1, 1], [1, 2]]         => 15
  [[2, 0], [0, 3]]         => 90

julia> endom = all_weights_endomorphisms_universal_bundle(Q,d,theta)
Dict{Vector{Vector{Int64}}, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [0, 5, -5, 0]
  [[2, 1], [0, 2]]         => [0, 10, -10, 0]
  [[1, 0], [1, 2], [0, 1]] => [0, 10, 15, -10, 0, 5, -15, -5, 0]
  [[1, 0], [1, 3]]         => [0, 15, -15, 0]
  [[1, 0], [1, 1], [0, 2]] => [0, 5, 10, -5, 0, 5, -10, -5, 0]
  [[1, 1], [1, 2]]         => [0, 5, -5, 0]
  [[2, 0], [0, 3]]         => [0, 5, -5, 0]

julia> all(maximum(endom[key]) < hn[key] for key in keys(hn))
true
```

The fact that all of these inequalities are satisfied allows to conclude that the higher cohomology of
``\mathcal{U}_i^\vee \otimes \mathcal{U}_j`` vanishes.

## Canonical decompositions and roots

QuiverTools provides a method to compute the canonical decomposition of a dimension vector.

Given a dimension vector ``d``, the canonical decomposition is a list of Shur roots ``\beta_i`` such that
``d = \sum_i \beta_i`` and ``\mathrm{ext}(\beta_i,\beta_j) = 0`` for all ``i \neq j``.

It is a theorem of Schofield that the canonical decomposition exists and is described by the condition above. It is a fact that if this is the canonical deconstruction of ``d``, then the general representation of dimension vector ``d`` decomposes as a direct sum of indecomposable representations of dimension vector ``\beta_i``.

```julia-repl
julia> canonical_decomposition(Q, d)
1-element Vector{Vector{Int64}}:
 [2, 3]

julia> canonical_decomposition(Q, [12,3])
6-element Vector{Vector{Int64}}:
 [1, 0]
 [1, 0]
 [1, 0]
 [3, 1]
 [3, 1]
 [3, 1]

julia> canonical_decomposition(Q, [12,4])
4-element Vector{Vector{Int64}}:
 [3, 1]
 [3, 1]
 [3, 1]
 [3, 1]
```

QuiverTools also implements computations of generic hom and ext for dimension vectors.

```julia-repl
julia> e = [1,2];

julia> generic_hom(Q, d, e)
0

julia> generic_hom(Q, e, d)
0

julia> generic_ext(Q, d, e)
4

julia> generic_ext(Q, e, d)
1
```

This allows to determine wether a root is real, imaginary isotropic or imaginary anisotropic.

```julia-repl
julia> ds = QuiverTools.all_subdimension_vectors([5,5])
6×6 Matrix{Vector{Int64}}:
 [0, 0]  [0, 1]  [0, 2]  [0, 3]  [0, 4]  [0, 5]
 [1, 0]  [1, 1]  [1, 2]  [1, 3]  [1, 4]  [1, 5]
 [2, 0]  [2, 1]  [2, 2]  [2, 3]  [2, 4]  [2, 5]
 [3, 0]  [3, 1]  [3, 2]  [3, 3]  [3, 4]  [3, 5]
 [4, 0]  [4, 1]  [4, 2]  [4, 3]  [4, 4]  [4, 5]
 [5, 0]  [5, 1]  [5, 2]  [5, 3]  [5, 4]  [5, 5]

julia> filter(d -> QuiverTools.is_real_root(Q, d), ds)
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [3, 1]
 [1, 3]

julia> filter(d -> QuiverTools.is_isotropic_root(Q, d), ds)
1-element Vector{Vector{Int64}}:
 [0, 0]

julia> filter(d -> QuiverTools.is_imaginary_root(Q, d), ds)
20-element Vector{Vector{Int64}}:
 [0, 0]
 [1, 1]
 [2, 1]
 [1, 2]
 [2, 2]
 [3, 2]
 [4, 2]
 [5, 2]
 [2, 3]
 [3, 3]
 [4, 3]
 [5, 3]
 [2, 4]
 [3, 4]
 [4, 4]
 [5, 4]
 [2, 5]
 [3, 5]
 [4, 5]
 [5, 5]
```

## Hodge polynomials

QuiverTools features an implementation of the Hodge polynomial of quiver moduli, if the base field is ``\mathbb{C}`` and the dimension vector is a coprime Schurian root.

```julia-repl
julia> Q = mKronecker_quiver(3); d = [2, 3]; theta = canonical_stability(Q, d);

julia> Hodge_polynomial(Q, d, theta)
x^6*y^6 + x^5*y^5 + 3*x^4*y^4 + 3*x^3*y^3 + 3*x^2*y^2 + x*y + 1

julia> Hodge_diamond(Q, d, theta)
7×7 Matrix{Int64}:
 1  0  0  0  0  0  0
 0  1  0  0  0  0  0
 0  0  3  0  0  0  0
 0  0  0  3  0  0  0
 0  0  0  0  3  0  0
 0  0  0  0  0  1  0
 0  0  0  0  0  0  1
```

Note that the ``i, j``-th entry of the matrix representing the Hodge diamond is ``h^{i,j}``. In other words, the point of the diamond is on the upper left side of the matrix.

This allows us to conclude that the Picard rank of the moduli space is 1.

```julia-repl
julia> Picard_rank(Q, d, theta)
1
```

For performance-oriented computations, one can use some theoretical results to get a slightly faster computation of the Hodge polynomial.

```julia-repl
julia> QuiverTools._Hodge_polynomial_fast(Q, d, theta)
q^6 + q^5 + 3*q^4 + 3*q^3 + 3*q^2 + q + 1
```

This skips some safety checks and returns the Hodge polynomial after the change of variables ``q = x*y``. This leads to about a 3% speedup in the computation...

```julia-repl
julia> using BenchmarkTools

julia> @benchmark Hodge_polynomial(Q, d, theta)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  117.500 μs … 109.571 ms  ┊ GC (min … max):  0.00% … 41.09%
 Time  (median):     123.188 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   230.718 μs ±   3.241 ms  ┊ GC (mean ± σ):  18.81% ±  1.36%

  	   ▂▄██▇▅▄▃▂▂                                                  
  ▂▄▇███████████▇▆▅▅▄▄▃▃▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  118 μs           Histogram: frequency by time          155 μs <

 Memory estimate: 136.06 KiB, allocs estimate: 3408.

julia> @benchmark QuiverTools._Hodge_polynomial_fast(Q, d, theta)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  114.625 μs … 110.717 ms  ┊ GC (min … max):  0.00% … 40.68%
 Time  (median):     120.042 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   219.409 μs ±   3.148 ms  ┊ GC (mean ± σ):  18.09% ±  1.28%

        ▄▅▆▆█▇▆▇▅▄▂▂▁▁                                           
  ▁▂▃▅▇███████████████▇█▆▆▆▄▅▄▄▃▃▃▃▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  115 μs           Histogram: frequency by time          139 μs <

 Memory estimate: 132.34 KiB, allocs estimate: 3327.
```


<!-- 
## Bundle library

QuiverTools contains (for now) a basic implementation of Bundle objects.

These are meant to be used as containers to perform computations with weights
or linearizations.

### Functionality

One can create a Bundle object by passing a list of weights and the rank to the constructor:

```julia-repl
julia> Bundle([1, 2, 3], 3)
Bundle of rank 3 with weights: [1, 2, 3]
```

`Note: the rank is not computed automatically from the weights. This is a design choice, to force a sanity check on the user.
Should there be the necessity to bypass this, the expected rank is the length of the list of weights.`

Direct sum, tensor product, wedge and symmetric product are implemented:

```julia-repl
julia> U = Bundle([1, 2, 3], 3); V = Bundle([4, 5], 2);
julia> U ⊕ V
Bundle of rank 5 with weights: [1, 2, 3, 4, 5]

julia> U ⊗ V
Bundle of rank 6 with weights: [5, 6, 6, 7, 7, 8]

julia> det(U)
Bundle of rank 1 with weights: [6]

julia> wedge(U, 0)
Bundle of rank 1 with weights: [0]
julia> wedge(U, 1)
Bundle of rank 3 with weights: [1, 2, 3]
julia> wedge(U, 2)
Bundle of rank 3 with weights: [3, 4, 5]
julia> wedge(U, 3)
Bundle of rank 1 with weights: [6]
julia> wedge(U, 4)
Bundle of rank 0 with weights: Int64[]

julia> symm(U,0)
Bundle of rank 1 with weights: [0]
julia> symm(U,1)
Bundle of rank 3 with weights: [1, 2, 3]
julia> symm(U,2)
Bundle of rank 6 with weights: [2, 3, 4, 4, 5, 6]
julia> symm(U,3)
Bundle of rank 10 with weights: [3, 4, 5, 5, 6, 7, 6, 7, 8, 9]
julia> symm(U,4)
Bundle of rank 15 with weights: [4, 5, 6, 6, 7, 8, 7, 8, 9, 10, 8, 9, 10, 11, 12]
```

A basic implementation of box products is also available:
  
```julia-repl

julia> U = Bundle([1, 2, 3], 3); V = Bundle([4, 5], 2);
julia> a = U ⊠ V
Bundle of rank 6 with weights: [[1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5]]
julia> b = V ⊠ V
Bundle of rank 4 with weights: [[4, 4], [4, 5], [5, 4], [5, 5]]



julia> a ⊗ b
Bundle of rank 24 with weights: [[5, 8], [5, 9], [6, 8], [6, 9], [5, 9], [5, 10], [6, 9], [6, 10], [6, 8], [6, 9], [7, 8], [7, 9], [6, 9], [6, 10], [7, 9], [7, 10], [7, 8], [7, 9], [8, 8], [8, 9], [7, 9], [7, 10], [8, 9], [8, 10]]

julia> a ⊕ b
Bundle of rank 10 with weights: [[1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5], [4, 4], [4, 5], [5, 4], [5, 5]]

julia> wedge(b,2)
Bundle of rank 6 with weights: [[8, 9], [9, 8], [9, 9], [9, 9], [9, 10], [10, 9]]

julia> wedge(b,3)
Bundle of rank 4 with weights: [[13, 13], [13, 14], [14, 13], [14, 14]]

julia> wedge(b,5)
Bundle of rank 0 with weights: Vector{Int64}[]

julia> symm(a,2)
Bundle of rank 21 with weights: [[2, 8], [2, 9], [3, 8], [3, 9], [4, 8], [4, 9], [2, 10], [3, 9], [3, 10], [4, 9], [4, 10], [4, 8], [4, 9], [5, 8], [5, 9], [4, 10], [5, 9], [5, 10], [6, 8], [6, 9], [6, 10]]
```

One use case is creating Bundle objects containing weights
of linearisations with respect to 1-PSs and then verifying Teleman inequality for
objects built with complex combinations of tensor, box and exterior products.

Another use case is defining line bundles on a projective space and then computing
the effects of direct summands, tensor products and exterior products on the
linearization. One such example is the computation of the Eagon-Northcott complex
for P^n, done below.

```julia-repl
julia> Q = mKronecker_quiver(2);
julia> U = [Bundle([-1], 1),Bundle([0], 1)]
2-element Vector{Bundle}:
  Bundle of rank 1, with weights: [-1]
  Bundle of rank 1, with weights: [0]

julia> EagonNorthcottcomplex(Q,U)
1-element Vector{Bundle}:
  Bundle of rank 1, with weights: [[-1, -1]]


julia> Q = mKronecker_quiver(5);
julia> EagonNorthcottcomplex(Q,U)
4-element Vector{Bundle}:
  Bundle of rank 10, with weights: [[-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1]]
  Bundle of rank 20, with weights: [[-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1]]
  Bundle of rank 15, with weights: [[-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1]]
  Bundle of rank 4, with weights: [[-1, -4], [-2, -3], [-3, -2], [-4, -1]]
```

Each Bundle object contains a rank and a list of weights.
 -->