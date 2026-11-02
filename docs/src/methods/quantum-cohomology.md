# Quasimap quantum cohomology

```@meta
CurrentModule = QuiverTools
DocTestSetup = :(using QuiverTools)
```

QuiverTools computes a finite polynomial presentation whose completion in the
effective Novikov monoid is the quasimap quantum cohomology of a smooth projective
quiver moduli space.
The construction uses [Zhang--Zhou, Theorem 2.22](https://arxiv.org/abs/2511.09875v2).
It works over the rational numbers, without enumerating quantum points or assuming
semisimplicity.

The input is an acyclic quiver with a nonzero, nonnegative, indivisible dimension
vector and an orthogonal, coprime stability parameter, with nonempty stable locus.
These hypotheses ensure that the effective group acts freely and has torsion-free
fundamental group.
Vertices of dimension zero are supported.
For a framed moduli object, apply the constructor to its [`total_space`](@ref).

## First example: projective space

The moduli space for the ``n``-Kronecker quiver and dimension vector ``(1,1)`` is
projective space of dimension ``n-1``.
Its hyperplane class satisfies ``h^n=q``.

```jldoctest
julia> A = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(3), [1, 1]));

julia> x = quantum_tautological_classes(A); q = only(novikov_variables(A));

julia> h = x[2][1] - x[1][1];

julia> iszero(quantum_normal_form(A, h^3 - q))
true

julia> gw_comparison(A)
:equal

julia> matrix, basis, ring = quantum_multiplication_matrix(A, h, [1]);

julia> matrix^3 == one(matrix)
true
```

Elements have Singular parents.
Use [`quantum_normal_form`](@ref) to reduce a polynomial expression in the presentation.
The matrix columns represent multiplication of the corresponding basis elements.
The basis is a standard monomial basis, not a specified Schubert basis.

## The classical limit and nilpotents

Setting every Novikov generator to zero gives the classical cohomology ring.
At this specialization the hyperplane is nilpotent.
This information is retained in the multiplication matrix:

```jldoctest
julia> A = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(3), [1, 1]));

julia> x = quantum_tautological_classes(A); h = x[2][1] - x[1][1];

julia> matrix, basis, ring = quantum_multiplication_matrix(A, h, [0]);

julia> length(basis), iszero(matrix^3), iszero(matrix^2)
(3, true, false)
```

Neither arrow weights nor the Weyl discriminant are inverted.
Consequently the presentation retains components with coincident Chern roots,
repeated eigenvalues, and nonreduced components.
For example, on ``\operatorname{Gr}(2,4)`` at ``q=1`` the hyperplane has a
two-dimensional zero eigenspace:

```jldoctest
julia> A = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(4), [1, 2]); chi=[1, 0]);

julia> c1, c2 = quantum_tautological_classes(A)[2]; q = only(novikov_variables(A));

julia> iszero(quantum_normal_form(A, c1^3 - 2c1*c2))
true

julia> iszero(quantum_normal_form(A, c1^4 - 3c1^2*c2 + c2^2 + q))
true

julia> matrix, basis, ring = quantum_multiplication_matrix(A, c1, [1]);

julia> length(basis), QuiverTools.Singular.Nemo.rank(matrix)
(6, 4)
```

## Several Novikov parameters

The effective degree monoid is computed integrally in
``\pi_1(G)=\mathbb Z^{\operatorname{supp}(\mathbf d)}/\mathbb Z\mathbf d``.
QuiverTools chooses its minimal additive generators and includes their binomial
relations, if any.
Thus Novikov generators are not silently treated as independent variables.
[`novikov_degrees`](@ref) returns representatives in the vertex lattice;
adding a multiple of the dimension vector does not change the degree.
Representatives are extended by zero at vertices outside the support.

Here is a presentation for ``\mathbb P^1\times\mathbb P^1``:

```jldoctest
julia> M = QuiverModuliSpace(three_vertex_quiver(2, 2, 0), [1, 1, 1]);

julia> A = quasimap_quantum_cohomology(M; chi=[1, 0, 0]);

julia> length(novikov_variables(A))
2

julia> x = quantum_tautological_classes(A); h = x[2][1] + x[3][1];

julia> matrix, basis, ring = quantum_multiplication_matrix(A, h, [1, 1]);

julia> length(basis), matrix^4 == 4matrix^2
(4, true)
```

The values supplied to specialization must respect all monoid relations.
All zeros and all ones always do so.
Exact integers and rationals are accepted; floating point values are rejected.
Some presentation degrees can describe base points without corresponding to
nonconstant curves in the moduli space.
For example, the one-arrow ``(1,1)`` moduli space is a point, but its quasimap
presentation still has a Novikov generator.
The degree group is the group of the GIT presentation, not an assumed copy of
the homology lattice of the quotient.

A four-vertex example has four minimal Novikov generators with a relation:

```jldoctest
julia> Q = Quiver([0 1 1 0; 0 0 1 1; 0 0 0 1; 0 0 0 0]);

julia> M = QuiverModuliSpace(Q, [1, 1, 1, 2], [6, 3, 1, -5]);

julia> A = quasimap_quantum_cohomology(M); q = novikov_variables(A);

julia> length(q), iszero(quantum_normal_form(A, q[1]*q[4] - q[2]*q[3]))
(4, true)
```

For this example, `[1, 1, 1, 2]` is not a valid specialization: it violates
``q_1q_4=q_2q_3`` and is rejected.

## A nonabelian example

For the ``3``-Kronecker moduli space with dimension vector ``(2,3)``, the classical
rank and the rank of the polynomial fiber at ``q=1`` are both ``13``:

```jldoctest
julia> A = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(3), [2, 3]));

julia> x = quantum_tautological_classes(A); h = -3x[1][1] + 2x[2][1];

julia> matrix, basis, ring = quantum_multiplication_matrix(A, h, [1]);

julia> length(basis), gw_comparison(A)
(13, :equal)
```

The normalization `chi` satisfies ``\sum_i\chi_i d_i=1`` and fixes the universal
bundles by imposing ``\sum_i\chi_i c_1(U_i)=0``.
It changes the Chern generators but not the underlying completed algebra.
The quantum constructor does not modify a moduli object's cached Chow ring or
the normalization selected for intersection theory.
The tests additionally compare this example with all five published quantum
relations of [Meng, Theorem 5.23](https://arxiv.org/abs/2412.15987).

## Comparison with Gromov--Witten quantum cohomology

The computed algebra is **quasimap** quantum cohomology.
[`gw_comparison`](@ref) returns `:equal` if every nonzero effective quasimap degree
has anticanonical degree at least two, and `:unknown` otherwise.
This is a sufficient test, not a characterization.
In the first case grading kills the correction
``\mu_\beta(z)=[zI_\beta(z)]_{z^{\geq0}}`` in
[Zhou's quasimap wall-crossing formula](https://arxiv.org/abs/1911.02745).
The three-relative-point construction agrees with three-marked ``0^+`` theory,
so its primary product agrees with the Gromov--Witten product.

The test runs on generators of the effective **quasimap** monoid.
Positivity merely on actual curves does not establish the required hypothesis.
For the canonical ``n``-Kronecker examples the determinant of the representation
is the ``n``th power of the primitive ample character.
Hence the test succeeds when ``n\geq2``.

Even when `gw_comparison(A) == :equal`, the generators represent **quantum
tautological classes**.
They need not be the classical Chern classes of the same names.
The comparison of products does not identify every polynomial expression with
its classical restriction.
No mirror transformation is computed when the comparison is `:unknown`.

## Polynomial presentations, completion, and specialization

The returned quotient is an exact **polynomial model**; its Novikov completion
is quasimap quantum cohomology.
It is not a calculation truncated at a chosen degree.
All specialization functions operate on this polynomial model.
They never attempt to evaluate an arbitrary formal power series at ``q=1``.

For general stability, the polynomial model can have components away from the
formal neighborhood of the Novikov origin.
A nonzero fiber can therefore differ from the completed theory and need not
have the classical rank.
No automatic interpretation as the ordinary quantum spectrum is made.
A positive grading, such as the canonical Kronecker grading, supplies the extra
control needed to identify the polynomial model with polynomial quantum
cohomology before specializing.

## Why finitely many cut degrees suffice

Expand vertex ``i`` into ``d_i`` torus vertices.
Every arrow from ``i`` to ``j`` contributes all weights ``\xi_{j,b}-\xi_{i,a}``,
with its original multiplicity.
A torus cocharacter is an integral potential ``k`` on these expanded vertices,
modulo common translation.
A generic section of its associated representation uses exactly the arrows
of nonnegative potential difference.
Thus effectivity is an exact stability test on that support.

For an effective potential, let ``S_l=\{v:k_v\leq l\}``.
Each nontrivial level cut ``-1_{S_l}`` is effective: its allowed-arrow support
contains the original allowed-arrow support.
Moreover,

```math
\begin{equation}
k-\max(k)\mathbf1=\sum_{l=\min(k)}^{\max(k)-1}-\mathbf1_{S_l}.
\end{equation}
```

All summands have compatible signs on every arrow weight.
Writing ``P_d`` and ``N_d`` for the positive and negative weight products, and
``s_d=(-1)^{\langle2\rho,d\rangle}``, their binomials consequently satisfy

```math
\begin{equation}
B_{d+e}=P_eB_d+s_dq^{\bar d}N_dB_e.
\end{equation}
```

It follows that the finitely many effective binary cuts generate all abelian
relations, not just those below a degree bound.
The implementation then removes redundant Novikov monoid generators and
computes the monoid kernel using auxiliary Laurent coordinates.
Only that kernel is transferred back; no Novikov variable is inverted in
the presentation.

Finally the Chern-root polynomial ring is free over the Weyl-invariant ring.
A coinvariant basis of size ``\prod_i d_i!`` therefore suffices for the factors
``g`` in the divided alternants ``\Delta^{-1}\operatorname{Alt}(gB_d)``.
Cuts in the same Weyl orbit give the same invariant ideal, since a Weyl
permutation of ``B_d`` can be transferred to the arbitrary multiplier ``g``
inside alternation, up to its sign.
It therefore suffices to retain one cut for each tuple of selected-vertex
counts in the original vertex blocks.

Each divided alternant is computed by composing adjacent divided differences
along a reduced word for the longest Weyl-group element.
For the discriminant convention above, an adjacent operator is
``f\mapsto(f-s_i f)/(\xi_{i+1}-\xi_i)``.
This uses ``\sum_i d_i(d_i-1)/2`` operators instead of summing over the full
Weyl group for each basis element.
Every polynomial division is exact, followed by descent to Chern generators.
No root or discriminant is inverted.

The cut enumeration is exponential in ``\sum_i d_i`` and alternation is
still factorial in the vertex dimensions through the coinvariant basis size.
The interface is general; large dimension vectors can still be impractical.
Use `verbose=true` to see the cut count, Weyl-orbit count, and relation-generation stage.

## Reference

```@docs
QuasimapQuantumCohomology
quasimap_quantum_cohomology
quasimap_presentation
novikov_variables
novikov_degrees
quantum_tautological_classes
quantum_normal_form
gw_comparison
specialize_quantum_parameters
quantum_multiplication_matrix
```
