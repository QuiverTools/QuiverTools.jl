# Quiver moduli

The main purpose of QuiverTools is to treat quiver moduli.
This package implements methods to study stable and semistable loci,
Harder--Narasimhan stratifications and quiver moduli spaces.

## Types

The following types are exported by QuiverTools:

```@docs
QuiverModuli
QuiverModuliSpace
QuiverModuliStack
```
We also have some basic constructors for well-known quiver moduli:

```@docs
kronecker_moduli
subspace_quiver_moduli
```

Black-box methods are provided to study some of their properties.

```@docs
is_nonempty
dimension
is_smooth
is_projective
index
motive
picard_rank
semistable_equals_stable
codimension_unstable_locus
```

## Stability

```@docs
has_semistables
has_stables
canonical_stability
is_coprime
slope
```

## Stratifications

QuiverTools implements methods to study the Harder--Narasimhan stratification
of the unstable locus of quiver representations, as well as the Luna stratification
of the properly semistable locus.

### Harder--Narasimhan types

```@docs
HNType
all_hn_types
is_hn_type
codimension_hn_stratum
is_amply_stable
```

### Luna types

A Luna type of a dimension vector ``\mathbf{d}`` for a stability parameter ``\theta`` is
an unordered sequence ``(\mathbf{d}^1, m_1), \dots, (\mathbf{d}^s, m_s)`` of dimension
vectors ``\mathbf{d}^k`` and positive multiplicities ``m_k`` with
``\sum_k m_k \mathbf{d}^k = \mathbf{d}``, all of the same ``\theta``-slope as
``\mathbf{d}``, and each ``\mathbf{d}^k`` admitting a ``\theta``-stable representation.
Luna types index the strata of the Luna stratification of the whole moduli space
``M^{ss}_\theta(Q, \mathbf{d})``: the open stratum is the stable locus (the trivial type
`Dict(d => [1])`), and the non-trivial types stratify the properly semistable locus. Each
stratum is described étale-locally by a local quiver built from the stable summands.

In QuiverTools a Luna type is represented by a [`LunaType`](@ref), which wraps a
dictionary whose keys are the distinct dimension vectors ``\mathbf{d}^k`` and whose values
are lists of the multiplicities with which they occur. For instance,
`Dict([1, 1] => [2, 1])` is the Luna type in which `[1, 1]` appears twice, once with
multiplicity `2` and once with multiplicity `1`. See [`LunaType`](@ref) for the precise
encoding.

```@docs
LunaType
all_luna_types
is_luna_type
dimension_of_luna_stratum
semisimple_moduli_space
```

## Hodge diamonds

```@docs
hodge_diamond
hodge_polynomial
betti_numbers
```
