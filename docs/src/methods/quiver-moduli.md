```@meta
DocTestSetup = quote
  using QuiverTools
end
```

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
codimension_singular_locus
fibre_dimension
is_flat
is_semismall
is_projective
index
motive
picard_rank
semistable_equals_stable
codimension_unstable_locus
```

## Stability

```@docs
has_semistable_representation
has_stable_representation
has_properly_semistable_representation
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
is_strongly_amply_stable
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
local_quiver_setting
local_structure
```

## Framed quiver moduli

A framed quiver moduli space ``M^{\Theta\text{-fr}}(Q, \mathbf{d}, \mathbf{n})`` comes
with a projection ``p`` to the quiver moduli space ``M^{\Theta}(Q, \mathbf{d})``.
The fibres of ``p`` over each Luna stratum of the base are described as the
nilpotent locus of a framed moduli space of the local quiver, following
[arXiv:2607.12895](https://arxiv.org/abs/2607.12895).

QuiverTools realizes it as a `FramedQuiverModuliSpace`,
which contains exactly the data of a `QuiverModuliSpace`,
together with a Boolean indicator of whether it is framed or coframed.
The nilpotent locus itself is wrapped in a `NilpotentLocus` structure.

From both of these new data types one can extract the underlying quiver moduli space,
respectively via `total_space(X::FramedQuiverModuliSpace)` and `ambient(X::NilpotentLocus)` -- see below.


```@docs
FramedQuiverModuliSpace
NilpotentLocus
```

```@docs
framing_vector
base
total_space
fibre
ambient
```

## Hodge diamonds

```@docs
hodge_diamond
hodge_polynomial
betti_numbers
```

## Intersection cohomology

If some proper subdimension vector has the slope of ``\mathbf{d}`` there are properly
semistable representations and the moduli space is singular, so the methods above no
longer apply. Its intersection cohomology is still computable, by the theorem of
[Meinhardt--Reineke](https://mathscinet.ams.org/mathscinet-getitem?mr=4000572)
identifying the Donaldson--Thomas invariants of the quiver with it. This needs the
stability parameter to be generic for the slope of ``\mathbf{d}``, and stable
representations to exist.

For instance, the 3-Kronecker quiver with dimension vector ``(2,2)`` has the ``(1,1)``'s
of the same slope, so its moduli space is singular. It turns out to have the intersection
cohomology of ``\mathbb{P}^5``:

```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 2]);

julia> intersection_poincare_polynomial(M)
L^5 + L^4 + L^3 + L^2 + L + 1
```

The quiver need not be acyclic. With one vertex, ``m`` loops and the trivial stability
parameter one gets the classical spaces of matrix invariants, which are affine, so that
the answer is then the Poincaré polynomial for cohomology with compact support.

```@docs
intersection_poincare_polynomial
intersection_betti_numbers
intersection_hodge_diamond
```

## Domokos reductions

[[Domokos](https://doi.org/10.4171/JCA/97)] introduces two operations on a quiver, a
dimension vector and a stability parameter which change the triple but preserve the
moduli space up to isomorphism: [`tau_reduction`](@ref) merges away a *large* vertex,
and [`sigma_reduction`](@ref) reflects a *small* source or sink. Iterating them
produces smaller moduli-isomorphic triples, and a pair from which no sequence of
reductions reaches anything smaller is *τσ-minimal*; this drives the finiteness of
quiver moduli spaces of a fixed dimension.

Stability parameters are expressed in QuiverTools' sign convention, which is opposite
to that of the paper; see the docstrings for the precise translation.

```@docs
is_large
is_small_source
is_small_sink
tau_reduction
sigma_reduction
tau_sigma_reduce
is_taus_minimal
```
