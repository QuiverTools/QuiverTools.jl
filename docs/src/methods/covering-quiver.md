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
free abelian group on ``Q_1``.

A *compatible dimension vector* for ``d \in \mathbb{N}^{Q_0}`` is a function
``\beta\colon Q_0 \times \mathbb{Z}^m \to \mathbb{N}`` with finite support such that
``\sum_\xi \beta(i, \xi) = d_i`` for every ``i \in Q_0``. The group ``\mathbb{Z}^m``
acts on compatible dimension vectors by
``s_\chi(\beta)(i, \xi) = \beta(i, \xi + \chi)``.

Every connected component of the natural-torus fixed locus of
``M^\theta(Q, d)`` is of the form
``F_\beta \cong M^{\hat\theta}(Q(w), \beta)``, where
``\hat\theta_{i, \xi} = \theta_i``, for a shift-equivalence class of compatible
``\beta``. A class contributes only when this lifted stable moduli space is
nonempty. This is the distinction between the candidates returned by
[`compatible_dimension_vectors`](@ref) and the actual components returned by
[`torus_fixed_components`](@ref); see Section 3 of
[[arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)],
building on Weist's localisation
(*Localization in quiver moduli spaces*, Represent. Theory 17 (2013), 382–425).

The dimension of the ``\chi``-weight space of the tangent space at the
fixed point is given by
[[Theorem 6.1, arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)]:

```math
\dim (T_{[M]} \mathcal M)_\chi
  = \delta_{\chi, 0} - \langle \beta, s_{-\chi}\beta\rangle_{Q(w)}.
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
nonzero_weights
```
