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

Black-box methods are provided to study some of their properties.

```@docs
is_nonempty
dimension
is_smooth
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
all_hn_types
is_hn_type
codimension_hn_stratum
is_amply_stable
```

### Luna types

```@docs
all_luna_types
is_luna_type
dimension_of_luna_stratum
```

## Hodge diamonds

```@docs
Hodge_diamond
Hodge_polynomial
Picard_rank
```
