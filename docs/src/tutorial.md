```@meta
DocTestSetup = quote
  using QuiverTools
end
```
# Tutorial

## Basic functionalities

To start using QuiverTools in the REPL, one first must import it.

```julia-repl
julia> using QuiverTools
```

Quivers can be built by passing the adjacency matrix to the `Quiver()` constructor:

```doctest
julia> Quiver([0 3; 0 0])
Quiver with adjacency matrix [0 3; 0 0]
```

The constructor accepts an optional string for naming the quiver:

```doctest
julia> MyQ = Quiver([0 3; 0 0], "My personal quiver")
My personal quiver, with adjacency matrix [0 3; 0 0]
```

QuiverTools has several constructors in place for many common examples. See
[Constructors](@ref) for all the special constructors.

```doctest
julia> kronecker_quiver(4)
4-Kronecker quiver, with adjacency matrix [0 4; 0 0]

julia> loop_quiver(5)
5-loop quiver, with adjacency matrix [5;;]

julia> subspace_quiver(3)
3-subspace quiver, with adjacency matrix [0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 0 0]

julia> three_vertex_quiver(1, 6, 7)
An acyclic 3-vertex quiver, with adjacency matrix [0 1 6; 0 0 7; 0 0 0]
```

Lastly, there is a handy constructor using strings:

```doctest
julia> Q = Quiver("1--2,1---3,2-4,2--3")
Quiver with adjacency matrix [0 2 3 0; 0 0 2 1; 0 0 0 0; 0 0 0 0]
```

Dimension vectors and stability parameters are represented by `AbstractVector{Int}`
objects, while under the hood these are encoded using the StaticArrays package.

```doctest
julia> Q = kronecker_quiver(3); d = [2,3];

julia> θ = canonical_stability(Q, d)
2-element StaticArraysCore.SVector{2, Int64} with indices SOneTo(2):
  9
 -6

julia> is_coprime(d, θ)
true
```

Here, `is_coprime()` checks if ``d`` is θ-coprime, i.e., if none of the
proper subdimension vectors ``0 \neq d' \nleq d`` satisfies ``\theta \cdot d' = 0``.

The bilinear Euler form relative to a quiver Q of any two vectors
in ``\mathbb{Z}^{Q_0}`` can be computed:

```doctest
julia> Q = kronecker_quiver(3); d = [2,2]; e = [3,4];

julia> euler_form(Q, d, e)
-10

julia> euler_form(Q, e, d)
-4
```

This allows to verify whether any given dimension vector belongs
to the fundamental domain of the problem.

The fundamental domain is the cone of dimension vectors in ``\mathbb{Z}^{Q_0}``
such that the symmetric Tits form is negative on all the simple roots, i.e.,
for all vertices i,

```math
(s_i, d) := \langle d, s_i\rangle + \langle s_i, d\rangle  \leq 0,
```

where ``s_i`` is the dimension vector with all entries set to ``0`` and the i-th
set to ``1``.

```doctest
julia> Q = kronecker_quiver(3); d = [2,2];

julia> in_fundamental_domain(Q, d)
true

julia> in_fundamental_domain(Q, [1,3])
false
```

## Stable and semistable dimension vectors

One can check if semistable, respectively stable representations
exist for a given dimension vector and stability parameter:

```doctest
julia> Q = kronecker_quiver(3); d = [2,3]; θ = [3,-2];

julia> has_semistable_representation(Q, d, θ)
true

julia> has_stable_representation(Q, d, θ)
true

julia> K2 = kronecker_quiver(2);

julia> has_stable_representation(K2, [2,2], [1,-1])
false

julia> has_semistable_representation(K2, [2,2], [1,-1])
true
```

One can also determine whether stable representations exist at all
for a given dimension vector by checking if it is a Schur root:

```doctest
julia> Q = kronecker_quiver(3); d = [2, 2];

julia> is_schur_root(Q, d)
true

julia> K2 = kronecker_quiver(2);

julia> is_schur_root(K2, d)
false
```

## Quiver Moduli

QuiverTools implements the abstract `QuiverModuli` type and the two concrete types
`QuiverModuliSpace` and `QuiverModuliStack`.

```doctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3])
Quiver moduli space defined as follows:
 - quiver: 3-Kronecker quiver
 - dimension vector: [2, 3]
 - stability parameter: [9, -6]
 - condition: semistable
```

Some special constructors are implemented for Kronecker moduli and subspace moduli.

```doctest
julia> kronecker_moduli(3, 3, 4)
Quiver moduli space defined as follows:
 - quiver: 3-Kronecker quiver
 - dimension vector: [3, 4]
 - stability parameter: [12, -9]
 - condition: semistable

julia> subspace_quiver_moduli(7, 5)
Quiver moduli space defined as follows:
 - quiver: 7-subspace quiver
 - dimension vector: [1, 1, 1, 1, 1, 1, 1, 5]
 - stability parameter: [5, 5, 5, 5, 5, 5, 5, -7]
 - condition: semistable
```


Several functionalities of QuiverTools are accessible either directly,
by passing a quiver, dimension vector, stability parameter etc,
or directly via these objects.
See the docstring of each method for more information and examples.

## Harder-Narasimhan types

QuiverTools provides methods to study the Harder-Narasimhan stratification
of the parameter space ``\mathrm{R}(Q,\mathbf{d})``.

```doctest
julia> Q = kronecker_quiver(3); M = QuiverModuliStack(Q, [2, 3], [3, -2]);

julia> all_hn_types(M)
8-element Vector{HNType}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> is_amply_stable(M)
true
```

The method `is_amply_stable()` determines whether
the codimension of the θ-semistable locus,

```math
\mathrm{R}^{\theta-sst}(Q,\mathbf{d})\subset\mathrm{R}(Q,\mathbf{d}),
```

is at least 2.

The method `all_hn_types()` provides a list of
all the Harder-Narasimhan types that appear in the problem.

The method `all_teleman_bounds()` computes the bounds
to apply Teleman quantization on the non-dense strata.
The output is a dictionary whose keys are the HN types
and whose values are the weights themselves.

```doctest
julia> Q = kronecker_quiver(3); M = QuiverModuliStack(Q, [2, 3], [3, -2]);

julia> all_teleman_bounds(M)
Dict{HNType{2}, Int64} with 7 entries:
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

```doctest
julia> Q = kronecker_quiver(3); M = QuiverModuliStack(Q, [2, 3]);

julia> hn = all_teleman_bounds(M)
Dict{HNType{2}, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 60
  [[2, 1], [0, 2]]         => 100
  [[1, 0], [1, 2], [0, 1]] => 100
  [[1, 0], [1, 3]]         => 360
  [[1, 0], [1, 1], [0, 2]] => 270
  [[1, 1], [1, 2]]         => 15
  [[2, 0], [0, 3]]         => 270

julia> endom = all_weights_endomorphisms_universal_bundle(M)
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

julia> does_teleman_inequality_hold(M)
true
```

The fact that all of these inequalities are satisfied allows
to conclude that the higher cohomology of
``\mathcal{U}_i^\vee \otimes \mathcal{U}_j`` vanishes.

## Canonical decompositions and roots

QuiverTools provides a method to compute
the canonical decomposition of a dimension vector.

Given a dimension vector ``d``, the canonical decomposition
is a list of Schur roots ``\beta_i`` such that
``d = \sum_i \beta_i`` and
``\mathrm{ext}(\beta_i,\beta_j) = 0`` for all ``i, j``.

It is a theorem of Schofield that the canonical decomposition exists and
is described by the condition above.
If this is the canonical decomposition of ``d``,
then the general representation of dimension vector ``d`` is
a direct sum of indecomposable representations of dimension vector ``\beta_i``.

```doctest
julia> Q = kronecker_quiver(3);

julia> canonical_decomposition(Q, [2, 3])
1-element Vector{Vector{Int64}}:
 [2, 3]

julia> canonical_decomposition(Q, [12, 3])
6-element Vector{Vector{Int64}}:
 [1, 0]
 [1, 0]
 [1, 0]
 [3, 1]
 [3, 1]
 [3, 1]

julia> canonical_decomposition(Q, [12, 4])
4-element Vector{Vector{Int64}}:
 [3, 1]
 [3, 1]
 [3, 1]
 [3, 1]
```

QuiverTools also implements computations of general hom and ext for dimension vectors.

```julia-repl
julia> e = [1, 2];

julia> general_hom(Q, d, e)
0

julia> general_hom(Q, e, d)
0

julia> general_ext(Q, d, e)
4

julia> general_ext(Q, e, d)
1
```

This allows to determine whether a root is real, imaginary isotropic
or imaginary anisotropic.

```doctest
julia> ds = QuiverTools.all_subdimension_vectors([5, 5])
6×6 Matrix{Vector{Int64}}:
 [0, 0]  [0, 1]  [0, 2]  [0, 3]  [0, 4]  [0, 5]
 [1, 0]  [1, 1]  [1, 2]  [1, 3]  [1, 4]  [1, 5]
 [2, 0]  [2, 1]  [2, 2]  [2, 3]  [2, 4]  [2, 5]
 [3, 0]  [3, 1]  [3, 2]  [3, 3]  [3, 4]  [3, 5]
 [4, 0]  [4, 1]  [4, 2]  [4, 3]  [4, 4]  [4, 5]
 [5, 0]  [5, 1]  [5, 2]  [5, 3]  [5, 4]  [5, 5]

julia> filter(d -> is_real_root(Q, d), ds)
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [3, 1]
 [1, 3]

julia> filter(d -> is_isotropic_root(Q, d), ds)
1-element Vector{Vector{Int64}}:
 [0, 0]

julia> filter(d -> is_imaginary_root(Q, d), ds)
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

QuiverTools features an implementation of the Hodge polynomial of quiver moduli,
if the base field is ``\mathbb{C}`` and the dimension vector is a coprime Schurian root.

```doctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> hodge_polynomial(M)
x^6*y^6 + x^5*y^5 + 3*x^4*y^4 + 3*x^3*y^3 + 3*x^2*y^2 + x*y + 1

julia> hodge_diamond(M)
7×7 Matrix{Int64}:
 1  0  0  0  0  0  0
 0  1  0  0  0  0  0
 0  0  3  0  0  0  0
 0  0  0  3  0  0  0
 0  0  0  0  3  0  0
 0  0  0  0  0  1  0
 0  0  0  0  0  0  1
```

Note that the ``i, j``-th entry of the matrix representing the Hodge diamond
is ``h^{i,j}``.
In other words, the point of the diamond is on the upper left side of the matrix.

This allows us to conclude that the Picard rank of the moduli space is 1.

```julia-repl
julia> picard_rank(M)
1
```

## Chow rings

QuiverTools allows to compute the Chow ring for a given quiver moduli space, as well as
the point class, the Todd class and the Euler characteristic of a vector bundle, given
its Chern character or its Chern class.

```doctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> CH = chow_ring(M; chi=[2, -1])
Singular polynomial quotient ring (QQ),(x1_1,x1_2,x2_1,x2_2,x2_3),(dp(5),C)

julia> QuiverTools.quotient_ideal(CH)
Singular ideal over Singular polynomial ring (QQ),(x1_1,x1_2,x2_1,x2_2,x2_3),(dp(5),C) with generators (2*x1_1 - x2_1, 33*x1_2*x2_3 - 7*x2_2*x2_3, 3*x1_2*x2_2 - x2_2^2 + x2_1*x2_3, x2_3^3, x2_2*x2_3^2, x2_1*x2_3^2, 9082*x2_2^2*x2_3 - 26539*x2_1*x2_3^2, 3*x2_1*x2_2*x2_3 - 22*x2_3^2, 3*x2_1^2*x2_3 - 8*x2_2*x2_3, 110*x2_2^3 - 153*x2_1*x2_2*x2_3 - 1188*x2_3^2, 5*x2_1*x2_2^2 - 5*x2_1^2*x2_3 - 144*x1_2*x2_3 + 6*x2_2*x2_3, x2_1^2*x2_2 - 48*x1_2^2 + 80*x1_2*x2_2 - 24*x2_2^2 + 18*x2_1*x2_3, x2_1^3 + 8*x1_2*x2_1 - 6*x2_1*x2_2 + 12*x2_3, x1_2*x2_1^2 - 12*x1_2^2 + 12*x1_2*x2_2 - 4*x2_2^2 + 4*x2_1*x2_3, 6*x1_2^2*x2_1 - 4*x1_2*x2_1*x2_2 + x2_1*x2_2^2 - x2_1^2*x2_3 - 12*x1_2*x2_3 + 2*x2_2*x2_3, 4*x1_2^3 - 4*x1_2*x2_1*x2_3 + x2_1*x2_2*x2_3 - 2*x2_3^2)

julia> u1, u2 = universal_bundle(M, 1), universal_bundle(M, 2)
(Bundle of rank 2, Bundle of rank 3)

julia> endom = dual(u1) * u2
Bundle of rank 6

julia> chern_class(u1)
x1_1 + x1_2 + 1

julia> chern_class(u2)
x2_1 + x2_2 + x2_3 + 1

julia> degree(u1)
57

julia> degree(u2)
3648

julia> ω = canonical_bundle(M)
Bundle of rank 1

julia> chern_class(ω)
-3//2*x2_1

julia> degree(ω)
41553

julia> euler_characteristic(ω)
1

julia> OO = structure_sheaf(M)
Bundle of rank 1

julia> euler_characteristic(OO)
1

julia> # As seen in `dual(u1) * u2` in the example above, tensor calculus is implemented to some extent:

julia> Λ = map(n -> exterior_power(u2, n), 0:4)
5-element Vector{Bundle}:
 Bundle of rank 1
 Bundle of rank 3
 Bundle of rank 3
 Bundle of rank 1
 Bundle of rank 0

julia> map(chern_character, Λ)
5-element Vector{Singular.spoly{Singular.n_Q}}:
 1
 -4//3*x1_2*x2_1 + 1//2*x2_1^2 + 1//2*x2_1*x2_2 - 5//36*x2_2^2 + 7//18*x2_1*x2_3 - 1//360*x2_2*x2_3 - 1//2160*x2_3^2 + x2_1 - x2_2 - 3//2*x2_3 + 3
 8*x1_2^2 - 8//3*x1_2*x2_1 + x2_1^2 + 3//2*x2_1*x2_2 - 29//36*x2_2^2 + 14//9*x2_1*x2_3 + 383//3960*x2_2*x2_3 + 23//432*x2_3^2 + 2*x2_1 - x2_2 - 9//2*x2_3 + 3
 8*x1_2^2 - 4//3*x1_2*x2_1 + 1//2*x2_1^2 + x2_1*x2_2 - 2//3*x2_2^2 + 5//3*x2_1*x2_3 + 76//165*x2_2*x2_3 + 76//135*x2_3^2 + x2_1 - 2*x2_3 + 1
 0

julia> ⨂ = map(n -> symmetric_power(u2, n), 0:3)
6-element Vector{Bundle}:
 Bundle of rank 1
 Bundle of rank 3
 Bundle of rank 6
 Bundle of rank 10
 Bundle of rank 15
 Bundle of rank 21

julia> map(chern_character, ⨂)
4-element Vector{Singular.spoly{Singular.n_Q}}:
 1
 -4//3*x1_2*x2_1 + 1//2*x2_1^2 + 1//2*x2_1*x2_2 - 5//36*x2_2^2 + 7//18*x2_1*x2_3 - 1//360*x2_2*x2_3 - 1//2160*x2_3^2 + x2_1 - x2_2 - 3//2*x2_3 + 3
 8*x1_2^2 - 40//3*x1_2*x2_1 + 3*x2_1^2 + 11//2*x2_1*x2_2 - 109//36*x2_2^2 + 70//9*x2_1*x2_3 + 31//3960*x2_2*x2_3 + 17//720*x2_3^2 + 4*x2_1 - 5*x2_2 - 33//2*x2_3 + 6
 88*x1_2^2 - 184//3*x1_2*x2_1 + 10*x2_1^2 + 53//2*x2_1*x2_2 - 269//12*x2_2^2 + 331//6*x2_1*x2_3 + 217//264*x2_2*x2_3 + 337//432*x2_3^2 + 10*x2_1 - 15*x2_2 - 157//2*x2_3 + 10
```

In fact, `Bundle` objects can also store the Teleman weights.
These are implemented for known bundles, behave well with respect to all
tensor calculus operations, and streamline the application of the Teleman quantization
theorem.

```doctest
julia> teleman_weights(u1)
Dict{HNType{2}, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [-15, -15]
  [[2, 1], [0, 2]]         => [-20, -20]
  [[1, 0], [1, 2], [0, 1]] => [-15, -25]
  [[1, 0], [1, 3]]         => [-45, -90]
  [[1, 0], [1, 1], [0, 2]] => [-45, -60]
  [[1, 1], [1, 2]]         => [0, -5]
  [[2, 0], [0, 3]]         => [-45, -45]

julia> teleman_weights(u2)
Dict{HNType{2}, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [-15, -15, -30]
  [[2, 1], [0, 2]]         => [-20, -30, -30]
  [[1, 0], [1, 2], [0, 1]] => [-25, -25, -30]
  [[1, 0], [1, 3]]         => [-90, -90, -90]
  [[1, 0], [1, 1], [0, 2]] => [-60, -75, -75]
  [[1, 1], [1, 2]]         => [0, -5, -5]
  [[2, 0], [0, 3]]         => [-60, -60, -60]

julia> teleman_weights(ω)
Dict{HNType{2}, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [90]
  [[2, 1], [0, 2]]         => [120]
  [[1, 0], [1, 2], [0, 1]] => [120]
  [[1, 0], [1, 3]]         => [405]
  [[1, 0], [1, 1], [0, 2]] => [315]
  [[1, 1], [1, 2]]         => [15]
  [[2, 0], [0, 3]]         => [270]

julia> teleman_weights(endom)
Dict{HNType{2}, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [0, 0, -15, 0, 0, -15]
  [[2, 1], [0, 2]]         => [0, -10, -10, 0, -10, -10]
  [[1, 0], [1, 2], [0, 1]] => [-10, -10, -15, 0, 0, -5]
  [[1, 0], [1, 3]]         => [-45, -45, -45, 0, 0, 0]
  [[1, 0], [1, 1], [0, 2]] => [-15, -30, -30, 0, -15, -15]
  [[1, 1], [1, 2]]         => [0, -5, -5, 5, 0, 0]
  [[2, 0], [0, 3]]         => [-15, -15, -15, -15, -15, -15]
```


The ample generator of the Picard group of `M`.
```doctest
julia> L = Bundle(M, chern_character_line_bundle(M, [3, -2]));

julia> map(i -> euler_characteristic(L^i), 0:5)
6-element Vector{Singular.n_Q}:
 1
 20
 148
 664
 2206
 5999
```

## Hochschild cohomology

For acyclic quivers, the Hochschild cohomologies are described in
[[Proposition 1.6, MR1035222](https://mathscinet.ams.org/mathscinet/relay-station?mr=1035222)].

QuiverTools provides a method to compute the only nontrivial one.

```doctest
julia> Q = kronecker_quiver(3);

julia> first_hochschild_cohomology(Q)
8
```

## Walls and chambers decomposition

QuiverTools provides an implementation of the walls and chambers decomposition
of the space of stability parameters in a quiver moduli problem,
as described in [[doi:10.1016/j.jpaa.2025.108153](https://doi.org/10.1016/j.jpaa.2025.108153)].

The provided methods are described in [Walls and Chambers decomposition](@ref).

### Technical note

The boolean function [`git_equivalent`](@ref) does not require an optional dependency.
Constructing the walls, chambers, and VGIT fan uses a *package extension*
that depends on the `Oscar` algebra system.

To use it, one must

- Ensure `Oscar` is installed in the current environment, by running `using Pkg; Pkg.add("Oscar")`; and
- Load the Oscar package, possibly without its interface, by running `import Oscar`.

If `import Oscar` is not run, `QuiverTools` will run it automatically,
but it may trigger recompilation and produce unwanted output.
For non-interactive workflows, it is best to run it at the beginning of the session.

### Displaying GIT fans

A stability parameter ``\theta`` for which there are semistable representations of
dimension vector ``d`` satisfies ``\theta \cdot d = 0``, so the GIT fan of a quiver with
``n`` vertices lies in the hyperplane ``d^\perp \subset \mathbb{R}^n`` of dimension ``n - 1``.
`Oscar.visualize` only draws objects in an ambient space of dimension at most 3.
The method [`project_to_d_perp`](@ref) rewrites the fan in coordinates of
``d^\perp \cong \mathbb{R}^{n-1}`` by forgetting one coordinate, which gives a picture in
the plane for quivers with 3 vertices and in space for quivers with 4 vertices.
This projection is a linear isomorphism, so it preserves the cones and how they meet,
but not lengths or angles.

For a quiver with 3 vertices, the GIT fan is drawn in the plane:

```julia-repl
julia> import Oscar

julia> Q = three_vertex_quiver(2, 3, 4); d = [1, 2, 2];

julia> F = project_to_d_perp(vgit_fan(Q, d), d);

julia> Oscar.ambient_dim(F), Oscar.n_maximal_cones(F)
(2, 3)

julia> Oscar.visualize(F)
```

For a quiver with 4 vertices, the GIT fan is drawn in space:

```julia-repl
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> F = project_to_d_perp(vgit_fan(Q, d), d);

julia> Oscar.ambient_dim(F), Oscar.n_maximal_cones(F)
(3, 4)

julia> Oscar.visualize(F)
```

The same works for the chambers and the walls, which are drawn together as a list:

```julia-repl
julia> Oscar.visualize(project_to_d_perp(vgit_chambers(Q, d), d))

julia> Oscar.visualize(project_to_d_perp(vgit_walls(Q, d), d))
```
