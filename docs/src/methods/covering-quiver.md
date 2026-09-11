# Covering quiver

For a finite quiver ``Q`` with arrow set ``Q_1 = \{a_1, \ldots, a_m\}``,
the *covering quiver* ``Q(w)`` is the (infinite) quiver

```math
\begin{aligned}
Q(w)_0 &= Q_0 \times \mathbb{Z}^m,\\
Q(w)_1 &= \bigl\{ (s(a_k), \xi) \to (t(a_k), \xi + e_k)
  \;\bigm|\; a_k \in Q_1,\ \xi \in \mathbb{Z}^m \bigr\},
\end{aligned}
```

where ``(e_k)_{k=1}^m`` is the standard basis of ``\mathbb{Z}^m``. The projection
``(i, \xi) \mapsto i`` is the universal abelian cover of ``Q`` with respect to the
free abelian group on ``Q_1``. Coordinates use the order returned by
[`arrows`](@ref), so a character vector and an arrow index always refer to
the same ordering.

A *compatible dimension vector* for ``d \in \mathbb{N}^{Q_0}`` is a function
``\beta\colon Q_0 \times \mathbb{Z}^m \to \mathbb{N}`` with finite support such that
``\sum_\xi \beta(i, \xi) = d_i`` for every ``i \in Q_0``. The group ``\mathbb{Z}^m``
acts on compatible dimension vectors by
``s_\chi(\beta)(i, \xi) = \beta(i, \xi + \chi)``.

Over an algebraically closed field, every connected component of the
natural-torus fixed locus of the stable moduli space ``M^\theta(Q, d)`` is of
the form
``F_\beta \cong M^{\hat\theta}(Q(w), \beta)``, where
``\hat\theta_{i, \xi} = \theta_i``, for a shift-equivalence class of compatible
``\beta``. A class contributes only when this lifted stable moduli space is
nonempty. This is the distinction between the candidates returned by
[`compatible_dimension_vectors`](@ref) and the actual components returned by
[`torus_fixed_components`](@ref); see
[[Theorem 3.1, Boos--Franzen](https://doi.org/10.1112/blms.12649)] and
[[Theorem 3.8, Weist](https://doi.org/10.1090/S1088-4165-2013-00436-3)].

For a semistable moduli space, `torus_fixed_components` requires the stable and
semistable loci to agree. Constructing `M` with `condition="stable"` explicitly
requests the fixed locus of the stable locus. Candidate enumeration is
combinatorial and can grow quickly with ``d`` and ``|Q_1|``; it is intended for
small and medium dimension vectors.

At a point of a stable fixed component, the dimension of the ``\chi``-weight
space of the tangent space is given by
[[Theorem 6.1, Boos--Franzen](https://doi.org/10.1112/blms.12649)]:

```math
\dim (T_{[M]} \mathcal M)_\chi
  = \delta_{\chi, 0} - \langle \beta, s_{-\chi}\beta\rangle_{Q(w)}.
```

Accordingly, [`weight_space_dimension`](@ref) and
[`tangent_weight_multiplicities`](@ref) accept the ambient moduli space and
verify that ``\beta`` has a nonempty stable lift. The latter includes the zero
character when the fixed component itself has positive dimension.

## Typical workflow

First enumerate candidates, then filter them using stability:

```julia
Q = kronecker_quiver(3)
M = QuiverModuliSpace(Q, [2, 3])

candidates = compatible_dimension_vectors(Q, M.d)
components = torus_fixed_components(M)

length(candidates)  # 55 connected-support shift classes
length(components)  # 13 nonempty stable fixed components

component = first(components)
tangent_weight_multiplicities(M, component.beta)
```

## Data structure

```@docs
CoveringDimVector
```

## Operations on covering dimension vectors

```@docs
shift_beta
covering_euler_form
extract_finite_subquiver
compatible_dimension_vectors
```

## Fixed components

```@docs
torus_fixed_components
```

## Tangent weight spaces at fixed points

```@docs
weight_space_dimension
tangent_weight_multiplicities
```

## References

- M. Boos and H. Franzen, *Weight spaces and attracting sets for torus actions
  on quiver moduli*, Bulletin of the London Mathematical Society **54** (2022),
  1658--1682. [doi:10.1112/blms.12649](https://doi.org/10.1112/blms.12649),
  [arXiv:2002.12049](https://arxiv.org/abs/2002.12049).
- A. D. King, *Moduli of representations of finite-dimensional algebras*,
  Quarterly Journal of Mathematics **45** (1994), 515--530.
  [doi:10.1093/qmath/45.4.515](https://doi.org/10.1093/qmath/45.4.515).
- T. Weist, *Localization in quiver moduli spaces*, Representation Theory
  **17** (2013), 382--425.
  [doi:10.1090/S1088-4165-2013-00436-3](https://doi.org/10.1090/S1088-4165-2013-00436-3),
  [arXiv:0903.5442](https://arxiv.org/abs/0903.5442).
