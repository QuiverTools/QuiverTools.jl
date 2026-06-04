# Chow rings

QuiverTools implements Chow rings computations for quivers and their representations.
There also is an implementation of `Bundle` objects in the Chow ring, which enable
tensor calculus in the Chow ring.

`Bundle` objects can also carry Teleman weights, and these behave well with respect to
tensor calculus.

A `QuiverModuliSpace` comes with an intermediary structure, `ChowRing`, which is
functionally just a container for various Chow ring data.

```@docs
QuiverTools.ChowRing
```

## Presentation of Chow rings

Chow rings must be initialized manually, passing a choice of linearization to the constructor.
If no linearization is passed, the default one is used.

```@docs
chow_ring
point_class
todd_class
chern_class_line_bundle
chern_character_line_bundle
integral
euler_characteristic
```

## Bundle objects

```@docs
Bundle
chern_character
chern_class
chern_classes
QuiverTools.degree
rank
teleman_weights
```

A line bundle descending from a linearization above can be computed.

```@docs
line_bundle
```

Some special bundles can be computed out of the box.

```@docs
QuiverTools.zero_sheaf
structure_sheaf
canonical_bundle
universal_bundle
```

## Tensor calculus

```@docs
dual
exterior_power
symmetric_power
det
```
