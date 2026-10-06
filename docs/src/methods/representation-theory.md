# Representation theory

QuiverTools offers various methods to study the representation theory of quivers,
their Schur roots, their canonical decomposition, and their Hochschild cohomology.

## Polynomial semi-invariant dimensions

`semi_invariant_dimension(Q, d, weight)` computes the dimension of a polynomial
semi-invariant weight space for an acyclic quiver.
The mathematical basis is the Schur-multiplicity description of these dimensions
in Derksen--Schofield--Weyman, Proposition 8
([MR2351613](https://mathscinet.ams.org/mathscinet-getitem?mr=2351613)).
The implementation expands those multiplicities by Cauchy's formula and the
Littlewood--Richardson rule; it does not construct semi-invariant polynomials.

Here is the weight convention and the coefficient formula used by the code.
Work over ``\mathbb C``, and write ``V_v=\mathbb C^{d_v}`` at vertex ``v``.
A coordinate function on an arrow ``a:i\to j`` transforms as
``V_i\otimes V_j^*``, so it has positive determinant weight at ``i`` and negative
weight at ``j``.
Choose a partition ``\lambda_a`` for each arrow, with at most
``\min(d_i,d_j)`` parts.
Cauchy's formula gives

```math
\begin{equation}
\mathbb C[\operatorname{Rep}(Q,\mathbf d)]
\cong
\bigoplus_{\boldsymbol\lambda}
\bigotimes_{a:i\to j}
\left(\mathbf S_{\lambda_a}(V_i)\otimes
\mathbf S_{\lambda_a}(V_j^*)\right).
\end{equation}
```

At vertex ``v``, let ``A_v(\mu)`` and ``B_v(\nu)`` be the Schur coefficients of
the outgoing and incoming products, respectively:

```math
\begin{equation}
\prod_{a:s(a)=v}s_{\lambda_a}
=\sum_\mu A_v(\mu)s_\mu,
\qquad
\prod_{a:t(a)=v}s_{\lambda_a}
=\sum_\nu B_v(\nu)s_\nu.
\end{equation}
```

Pad partitions with zeroes to length ``d_v``.
The multiplicity of ``\det^{\sigma_v}`` at ``v`` is

```math
\begin{equation}
m_v(\boldsymbol\lambda,\sigma_v)
=
\sum_{\substack{\mu,\nu\\
\mu_k-\nu_k=\sigma_v\ \text{for }1\leq k\leq d_v}}
A_v(\mu)B_v(\nu).
\end{equation}
```

The returned dimension is the sum of the products of these local
multiplicities over all arrow partitions:

```math
\begin{equation}
\dim \mathrm{SI}(Q,\mathbf d)_\sigma
=\sum_{\boldsymbol\lambda}\prod_{v\in Q_0}
m_v(\boldsymbol\lambda,\sigma_v).
\end{equation}
```

In particular, a contributing assignment must satisfy
``\sum_{a:s(a)=v}|\lambda_a|-\sum_{a:t(a)=v}|\lambda_a|=\sigma_v d_v``
at each vertex.
For the three-arrow Kronecker quiver with ``\mathbf d=(1,1)`` and
``\sigma=(2,-2)``, an arrow partition is just a nonnegative integer.
The three integers must sum to ``2``, and every such choice has local
multiplicity one.
There are ``\binom{4}{2}=6`` choices, so the function returns ``6``.
The search assigns arrow partitions in topological order and uses this
equality to restrict their sizes.
At a pure source or sink, the required Schur coefficient is rectangular.
The implementation then uses rectangle complements; it also uses the
rectangular-square identity and Pieri's rule when the factor shapes allow it.
These simplifications depend on partitions rather than on a quiver family.
For a partition ``\lambda\subseteq(w^r)``, define its rectangle complement by
``\lambda_i^c=w-\lambda_{r+1-i}``.
The two-factor coefficient of the rectangle is

```math
\begin{equation}
[s_{(w^r)}](s_\lambda s_\mu)
=
\begin{cases}
1,&\mu=\lambda^c,\\
0,&\text{otherwise}.
\end{cases}
\end{equation}
```

For an equal rectangular pair, the LR rule gives the multiplicity-free identity

```math
\begin{equation}
s_{(w^h)}^2
=
\sum_{\tau\subseteq(w^h)}
s_{(w+\tau_1,\ldots,w+\tau_h,w-\tau_h,\ldots,w-\tau_1)}.
\end{equation}
```

Terms with too many parts for the vertex rank are omitted.
When a factor is a row, Pieri's rule gives

```math
\begin{equation}
s_\lambda s_{(k)}
=\sum_{\mu/\lambda\text{ a horizontal }k\text{-strip}}s_\mu.
\end{equation}
```

The number of possible arrow partitions can still make a calculation expensive.

The alternative `method=:reciprocity` uses Derksen--Weyman, Theorem 1 and
Corollary 1
([MR1758750](https://mathscinet.ams.org/mathscinet-getitem?mr=1758750)).
If `E=euler_matrix(Q)` and ``\sigma=E^{\mathsf T}\alpha`` for a nonnegative
dimension vector ``\alpha``, then, provided ``\sigma\cdot\mathbf d=0``,

```math
\begin{equation}
\dim \mathrm{SI}(Q,\mathbf d)_\sigma
=\dim \mathrm{SI}(Q,\alpha)_{-E\mathbf d}.
\end{equation}
```

The code solves for ``\alpha`` in topological order and evaluates the right-hand
side with the same Schur calculation.
The spanning theorem implies that a nonzero weight space requires a nonnegative
``\alpha``, so the method returns zero when this solution has a negative entry.
Reciprocity can be slower when ``\alpha`` is larger than ``\mathbf d``.
For an amply stable quiver moduli space, the semi-invariant weight space also
computes global sections of its descended line bundle.
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
