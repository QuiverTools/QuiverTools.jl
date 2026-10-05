# Representation theory

QuiverTools offers various methods to study the representation theory of quivers,
their Schur roots, their canonical decomposition, and their Hochschild cohomology.

`semi_invariant_dimension(Q, d, weight)` computes any polynomial semi-invariant
weight-space dimension for an acyclic quiver.
It uses Cauchy decomposition on the arrows and Littlewood--Richardson coefficients
at the vertices.
At a pure source or sink, the determinant weight selects a rectangular
Schur coefficient. The calculation splits its factors and pairs complementary
partitions; rectangular Schur squares and Pieri's rule speed up the resulting
products whenever those shapes occur. These rules depend on the local
partitions, not on a particular quiver family. The general partition search
can still be expensive when many arrow partitions remain possible.
The optional `method=:reciprocity` evaluates the dual weight space from
Derksen--Weyman reciprocity (MR1758751, Corollary 1).
It is useful when the reciprocal dimension vector gives a smaller calculation.
For an amply stable quiver moduli space, the corresponding semi-invariant space
also computes global sections of its descended line bundle.
In particular, `weight=canonical_stability(Q, d)` gives anticanonical sections
under that hypothesis.

```@docs
euler_form
euler_matrix
is_schur_root
is_real_root
is_imaginary_root
is_isotropic_root
general_ext
general_hom
all_general_subdimension_vectors
is_general_subdimension_vector
canonical_decomposition
in_fundamental_domain
first_hochschild_cohomology
semi_invariant_dimension
```
