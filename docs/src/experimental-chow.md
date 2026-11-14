# Experimental parabolic Chow relations

!!! warning "Experimental"
    This is an LLM-suggested optimization.
    The proof and code need careful checking.

The default algorithm is `:weyl`.
Use the parabolic algorithm on a fresh moduli object:

```jldoctest
julia> using QuiverTools

julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> CH = chow_ring(M; algorithm=:experimental_parabolic);

julia> numbers = chern_numbers(M); (numbers["c_1^6"], numbers["c_6"])
(41553, 13)

julia> chow_ring(M; algorithm=:weyl) === CH
true
```

An existing cached ring is reused when the normalization `chi` is unchanged,
regardless of `algorithm` or `adaptive`.
Changing `chi` rebuilds the ring and clears the cached point and Todd classes.
To compare presentations, construct a fresh moduli object for each algorithm.
Bundle and integration methods use whichever presentation is stored on that object.

Use `adaptive=false` to compute the full parabolic generating set in the relevant degrees,
without early stopping.
Both experimental modes impose the dimension truncation proved below.
The existing stability checks are retained;
`unsafe=true` makes the caller responsible for the geometric hypotheses.

Assume a nonempty, smooth, projective, fine moduli space of stable representations,
with stable equal to semistable and the quiver acyclic on the support of $d$.
Theta-coprimality is a sufficient condition.
Zero-dimensional vertex spaces are removed for the independent Betti computation.
The reduction applies in any dimension and for any stability parameter satisfying these assumptions.

## Choosing an algorithm

| Keywords | Relations | Completeness check |
|:--|:--|:--|
| `algorithm=:weyl` | Full Weyl-module basis, followed by elimination | Existing tautological presentation |
| `algorithm=:experimental_parabolic` | Parabolic Schur bases, in increasing Schur degree | Independent Betti numbers |
| `algorithm=:experimental_parabolic, adaptive=false` | Full parabolic Schur bases in the relevant degrees | Exhaustion of those bases |

There is no automatic selection.
The reduction in the number of basis elements helps when the vertex ranks are large.
It does not by itself predict the running time:
the dimension truncation can introduce many monomials, and adaptive mode also computes Betti numbers.
For example, the Weyl algorithm is faster on the 3-Kronecker sixfold above,
while the parabolic algorithm is faster on the measured rank $(5,5,2)$ sevenfolds.
The benchmark report below gives both improvements and regressions.

The updated measurements also test the proposed score $S=FW^2/B$,
where $W=\prod_i d_i!$, $F$ counts minimal forbidden vectors,
and $B$ counts the minimal weighted monomials above the dimension before normalization.
A single cutoff at 32 does not predict the winner reliably near the boundary.
Five of six additional examples expose warm parabolic wins below 32,
while a cache-cleared warm parabolic call on 5-Kronecker $(2,5)$ exceeds 30 seconds at score 19.34.
Weyl takes about 10.45 seconds there.
The current data suggest Weyl below 16 and parabolic above 32 as provisional starting points,
with scores 16–32 left for direct comparison.
These are timing observations for the measured implementations, rather than a general complexity bound.
First-call choices can differ because of compilation.
The full tables, the one-call warm bound, and the preselected checks are in the benchmark report.

## Implementation

The implementation follows these steps:

1. Keep the existing minimal forbidden vectors and normalization of the universal bundles.
2. Add the minimal monomial generators in degrees above the dimension.
3. For each forbidden vector, enumerate the rectangular partition bases at its vertices.
   Their Cartesian product is lazy.
4. Evaluate the parabolic and remaining divided differences.
   Coefficients stay in Singular's rational field throughout.
5. Normalize nonzero relations, discard duplicates using polynomial equality,
   and express them in elementary symmetric polynomials.
6. Compute a Gröbner basis.
   In adaptive mode, count standard monomials and compare with the independently computed Betti numbers;
   continue to the next Schur degree if they differ.

The internal routines are documented in the API reference:
[`QuiverTools.__chow_divided_difference`](@ref),
[`QuiverTools.__chow_divided_alternant`](@ref),
[`QuiverTools.__chow_parabolic_monomials`](@ref),
[`QuiverTools.__chow_parabolic_alternant`](@ref),
[`QuiverTools.__chow_elementary_expression`](@ref),
[`QuiverTools.__chow_truncation_relations`](@ref), and
[`QuiverTools.__chow_graded_dimensions`](@ref).
[`QuiverTools.__chow_ring_parabolic`](@ref) assembles the presentation and returns the stopping data used by the tests.
The public quiver-level `chow_ring` method returns its usual presentation tuple;
the moduli-level method stores the presentation on the moduli object.

## Parabolic reduction

Work over the rationals, and write

```math
\begin{equation}\label{equation:root-ring}
R=\mathbb Q[\xi_{i,r}],\qquad
W=\prod_i S_{d_i},\qquad A=R^W.
\end{equation}
```

For a forbidden dimension vector $e$, put $a_i=e_i$, $b_i=d_i-e_i$, and

```math
\begin{equation}\label{equation:forbidden-polynomial}
H=\prod_i(S_{a_i}\times S_{b_i}),\qquad
f_e=\prod_{\alpha:i\to j}\prod_{r=1}^{a_i}\prod_{s=a_j+1}^{d_j}
(\xi_{j,s}-\xi_{i,r}).
\end{equation}
```

In particular, $f_e$ is $H$-invariant, including with parallel arrows.
For $G=W,H$, let $\Delta_G$ be the product of
$\xi_{i,s}-\xi_{i,r}$ for $r<s$ within its blocks, and set

```math
\begin{equation}\label{equation:alternating-transfer}
p_G(g)=\frac{\mathrm{Alt}_G(g)}{\Delta_G},\qquad
\mathrm{Alt}_G(g)=\sum_{w\in G}\mathrm{sgn}(w)w(g).
\end{equation}
```

These are polynomial-valued operators: an alternating polynomial is divisible by
each of the corresponding differences, hence by their product.
Define the relative transfer, using representatives for the cosets $wH$ in $W/H$, by

```math
\begin{equation}\label{equation:relative-transfer}
T(h)=\frac{1}{\Delta_W}
\sum_{w\in W/H}\mathrm{sgn}(w)w(\Delta_Hh),\qquad h\in R^H.
\end{equation}
```

**Claim.** The $A$-module $p_W(f_eR)$ is generated by
$T(f_es_\lambda)$, where

```math
\begin{equation}\label{equation:schur-module}
s_\lambda=\prod_i s_{\lambda_i}(\xi_{i,1},\ldots,\xi_{i,a_i}),\qquad
\lambda_i\subseteq(b_i^{a_i}).
\end{equation}
```

There are $\prod_i\binom{d_i}{e_i}$ such generators, rather than
$\prod_i d_i!$ generators coming from a basis of $R$ over $A$.

**Proof.** Decomposing the alternating sum into cosets $wH$, and using the
$H$-invariance of $f_e$, gives, for every $b\in R$,

```math
\begin{equation}\label{equation:transfer-factorization}
p_W(f_eb)=T\bigl(f_ep_H(b)\bigr).
\end{equation}
```

The map $p_H:R\to R^H$ is surjective, since for every $h\in R^H$,

```math
\begin{equation}\label{equation:surjective-parabolic-transfer}
p_H\left(\frac{\Delta_Hh}{|H|}\right)=h.
\end{equation}
```

For a single block of rank $a+b$, the invariant ring
$\mathbb Q[x_1,\ldots,x_{a+b}]^{S_a\times S_b}$ is a free module over the
fully symmetric ring, with basis $s_\lambda(x_1,\ldots,x_a)$ for
$\lambda\subseteq(b^a)$.
The invariant ring has the presentation by the Chern classes of two bundles of
ranks $a,b$, with the identity $c(U)c(V)=c(E)$ over the ring of the Chern classes of $E$.
This is the equivariant Chow presentation of the Grassmann bundle of $a$-planes
in a universal rank $a+b$ bundle.
The Grassmann bundle formula gives a free module over the base, with the Schur
classes of the rank $a$ tautological bundle indexed by that rectangle as a basis
(up to the harmless dual-bundle signs).
Applying the formula to finite-dimensional Grassmannian approximations of the
classifying spaces gives the invariant-ring statement in every degree.
Taking the product over vertices gives the displayed $A$-basis of $R^H$.
Both $p_W$ and $T$ are $A$-linear.
Surjectivity of $p_H$ and the factorization therefore prove the equality of modules
and hence of the ideals they generate. $\square$

Franzen's tautological presentation identifies the Chow ideal with the sum of
these modules, together with the linear relation normalizing the universal bundles.
The existing selection of minimal forbidden vectors is unchanged.
The freedom to choose bases adapted to the forbidden vector already occurs in
[Franzen, Theorem 3.1 and its following remark, arXiv:1307.3066](https://arxiv.org/abs/1307.3066).
See also [Belmans–Franzen, Section 2.2, arXiv:2307.01711](https://arxiv.org/abs/2307.01711).

## Realization by divided differences

For each vertex use monomials with exponents

```math
\begin{equation}\label{equation:parabolic-monomials}
(\lambda_{i,1}+a_i-1,\ldots,\lambda_{i,a_i},
 b_i-1,\ldots,0).
\end{equation}
```

The determinantal definition of Schur polynomials shows that $p_H$ sends these
to the basis above, times the common unit sign
$(-1)^{\sum_i(\binom{a_i}{2}+\binom{b_i}{2})}$.
Thus it suffices to evaluate $p_W(f_em_\lambda)$.
The implementation first computes $p_H(m_\lambda)$ and multiplies by $f_e$,
then applies the remaining divided differences.

With the positive differences used above, the simple operator is

```math
\begin{equation}\label{equation:simple-divided-difference}
\partial_k(g)=\frac{g-s_k(g)}{x_{k+1}-x_k}.
\end{equation}
```

The longest-word operator is $p_W$.
Length additivity in $w_0(W)=(w_0(W)w_0(H))w_0(H)$ factors it through $p_H$.
In a vertex block split as $(a,b)$, the code constructs the remaining reduced word
by concatenating $(b+j-1,b+j-2,\ldots,j)$ for $j=1,\ldots,a$,
and applies operators in reverse word order.
Its length is $ab$, the length difference of the two longest elements.
The tests compare this factorization to the full operator for every split of rank five,
and compare the full operator to an explicit permutation sum.

The symmetric outputs are converted to elementary symmetric polynomials by
successively canceling the lexicographically leading monomial.
For a symmetric polynomial, its leading exponents in each block are weakly decreasing.
The product of elementary symmetric polynomials with exponents given by successive
differences has exactly that leading monomial, with coefficient one.
Subtraction strictly lowers the leading monomial in a fixed finite degree, proving
termination and the claimed expression without elimination.

## Truncation and stopping

Let $N=1-\langle d,d\rangle$ be the dimension in the stated geometric setting,
and let $K\subset A$ be the true Chow ideal, including the normalization relation.
Every homogeneous class above degree $N$ vanishes.
Therefore adjoining all monomials of weighted degree greater than $N$ leaves the
Chow quotient unchanged.
The code adds just their minimal monomial generators;
these have degrees between $N+1$ and $N+\max_i d_i$.
A larger minimal generator would still have degree greater than $N$ after removing
any variable dividing it, contradicting minimality.

For a fixed forbidden vector, write $F_e=\deg(f_e)$,
$D_H=\deg(\Delta_H)$ and $D_W=\deg(\Delta_W)$.
The monomial lifting a Schur basis element of degree $s$ has degree $D_H+s$.
Its divided alternant therefore has degree $F_e+D_H+s-D_W$, if it is nonzero.
The rectangle bases have maximum Schur degree $\sum_i e_i(d_i-e_i)$,
so only the following degrees can contribute:

```math
\begin{equation}\label{equation:schur-degree-bound}
0\leq s\leq
\min\left\{\sum_i e_i(d_i-e_i),\ N+D_W-F_e-D_H\right\}.
\end{equation}
```

A negative upper bound discards that forbidden vector entirely.
This bound applies to the Schur degree, rather than the degree of its lifting monomial,
and can exceed $N$ because the transfer lowers degree.
It is used in both experimental modes.

At an intermediate Schur degree, every computed relation is already proved to
lie in $K$.
Thus its ideal $J$ satisfies $J\subseteq K$, and there is a graded surjection

```math
\begin{equation}\label{equation:certified-surjection}
A/J\twoheadrightarrow A/K=\mathrm{CH}^*(M)_{\mathbb Q}.
\end{equation}
```

For these smooth projective fine acyclic moduli, the rational cycle map identifies
the Chow groups with even cohomology.
The existing Reineke computation of Betti numbers is independent of the relations.
The code counts standard monomials of a Gröbner basis in each weighted degree
$0,\ldots,N$ and stops only when their counts equal those even Betti numbers.
The surjection is then an isomorphism in each degree;
both sides vanish above $N$ by construction.
This proves that stopping early recovers the entire ideal.
Equality of Hilbert functions without the established inclusion $J\subseteq K$
would not suffice.
The implementation raises an error if any count is too small or the full relevant
set of relations does not attain the expected counts.

## Benchmarks

The timings include the parabolic reduction, divided differences, conversion to
elementary symmetric polynomials, truncation, and the Betti computation.
The smaller reference examples are slower with this algorithm.
See [the measurements and reproduction instructions](https://github.com/QuiverTools/QuiverTools.jl/blob/experimental/llm-parabolic-chow/benchmark/parabolic-chow.md).

The Hilbert check proves completeness once the relations are known to be valid.
Their validity depends on the transfer calculation above and its implementation;
the operator and whole-ideal tests check this on the listed examples.
