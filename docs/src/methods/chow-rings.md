# Chow rings

QuiverTools now encodes quiver-moduli intersection theory directly in Oscar's
experimental IntersectionTheory backend.

The quiver-specific work done by QuiverTools is:

- building the tautological Chow-ring presentation,
- choosing and caching the linearization,
- attaching the tautological bundles and tangent bundle,
- attaching Teleman weights to the standard quiver bundles.

The generic intersection-theory API comes from Oscar. In particular,
`chow_ring(M)` returns an `Oscar.AbstractVariety`, and then one uses Oscar's
`chow_ring`, `point_class`, `todd_class`, `integral`, `tautological_bundles`,
and `tangent_bundle` functions on that object.

```julia-repl
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> X = chow_ring(M; chi=[-1, 1])
AbstractVariety of dim 6

julia> Oscar.point_class(X)
x23^2
```

## Variety construction

```@docs
chow_ring
```

## Quiver bundles

```@docs
line_bundle
structure_sheaf
canonical_bundle
universal_bundle
chern_classes
teleman_weights
QuiverTools.degree
```

## Tensor calculus

For tensor operations such as `dual`, `det`, `exterior_power`, and
`symmetric_power`, use Oscar's bundle algebra directly on the returned
`AbstractBundle` objects.
