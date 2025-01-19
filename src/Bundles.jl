###############################################################################
# Tensor calculus on quiver moduli
###############################################################################

import Base: *, +, -, ^

export chern_character, chern_class, chern_classes, dual, exterior_power, symmetric_power,
  det, canonical_bundle, universal_bundle, degree, rank

_has_chern_data(F::Bundle) = isdefined(F, :chern_character) || isdefined(F, :chern_class)

function chern_character(F::Bundle)
  !_has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_character) &&
    setfield!(F, :chern_character, _chern_character_from_classes(F))
  return F.chern_character
end

function chern_classes(F::Bundle)
  !_has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, _chern_classes_from_character(F))
  return F.chern_class
end

function chern_class(F::Bundle)
  !_has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, _chern_classes_from_character(F))
  return sum(values(chern_classes(F)))
end

function chern_class(F::Bundle, k)
  !_has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, _chern_classes_from_character(F))
  return chern_classes(F)[k]
end

function teleman_weights(F::Bundle)
  !isdefined(F, :teleman_weights) && throw(ArgumentError("Bundle has no weights."))
  return F.teleman_weights
end

rank(F::Bundle) = F.rank
chow_ring(F::Bundle) = F.parent.ring
variety(F::Bundle) = F.parent.parent
structure_sheaf(M::QuiverModuliSpace) = Bundle(M, 1)

##############################
# Operations on Bundle objects
##############################

"""
    dual(F::Bundle)

Return the dual bundle of `F`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M))
Bundle of rank 1

julia> chern_character(dual(td))
-x21 + 1
```

On our favourite 6-fold:

```jldoctest
julia> Q, d = kronecker_quiver(3), [2, 3];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=[-1, 1]);

julia> td = Bundle(M, todd_class(M));

julia> chern_character(td)
-17//8*x12*x21 + x21^2 + 823//360*x12*x22 - 823//1080*x22^2 + 553//1080*x21*x23 - 77//60*x22*x23 + x23^2 + 5//12*x12 - 3//2*x21 + 9//8*x23 + 1


julia> tdd = dual(td)
Bundle of rank 1

julia> chern_character(tdd)
17//8*x12*x21 + x21^2 + 823//360*x12*x22 - 823//1080*x22^2 + 553//1080*x21*x23 + 77//60*x22*x23 + x23^2 + 5//12*x12 + 3//2*x21 - 9//8*x23 + 1
```
"""
function dual(F::Bundle)
  return Bundle(F.parent, adams(F, -1))
end

*(n::Int, F::Bundle) = Bundle(F.parent, n * chern_character(F))
*(F::Bundle, n::Int) = n * F
^(F::Bundle, n::Int) = Bundle(F.parent, chern_character(F)^n)

# direct sum, quotient and tensor product
+(F::Bundle, G::Bundle) =
  if F.parent == G.parent
    Bundle(F.parent, chern_character(F) + chern_character(G))
  else
    throw(DomainError("Different Chow rings."))
  end
-(F::Bundle, G::Bundle) =
  if F.parent == G.parent
    Bundle(F.parent, chern_character(F) - chern_character(G))
  else
    throw(DomainError("Different Chow rings."))
  end
*(F::Bundle, G::Bundle) =
  if F.parent == G.parent
    Bundle(F.parent, chern_character(F) * chern_character(G))
  else
    throw(DomainError("Different Chow rings."))
  end

"""
    exterior_power(F::Bundle, k::Int)

Return the `k`-th exterior power of `F`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M))
Bundle of rank 1

julia> F = 3*td
Bundle of rank 3

julia> chern_character(F)
3*x21 + 3

julia> W = map(i -> exterior_power(F, i), 0:4);

julia> map(w -> (rank(w), chern_character(w)), W)
5-element Vector{Tuple{Int64, Singular.spoly{Singular.n_Q}}}:
 (1, 1)
 (3, 3*x21 + 3)
 (3, 6*x21 + 3)
 (1, 3*x21 + 1)
 (0, 0)
```
"""
function exterior_power(F::Bundle, k::Int)
  return Bundle(F.parent, _chern_characters_wedge(F, k)[end])
end
det(F::Bundle) = exterior_power(F, rank(F))

"""
    symmetric_power(F::Bundle, k::Int)

Return the `k`-th symmetric power of `F`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M));

julia> F = 3*td
Bundle of rank 3

julia> chern_character(F)
3*x21 + 3

julia> W = map(i -> symmetric_power(F, i), 0:4);

julia> map(w -> (rank(w), chern_character(w)), W)
5-element Vector{Tuple{Int64, Singular.spoly{Singular.n_Q}}}:
 (1, 1)
 (3, 3*x21 + 3)
 (6, 12*x21 + 6)
 (10, 30*x21 + 10)
 (15, 60*x21 + 15)
```
"""
function symmetric_power(F::Bundle, k::Int)
  return Bundle(F.parent, _chern_characters_symmetric(F, k)[end])
end

function homogeneous_components(M::QuiverModuliSpace, x)
  n = dimension(M)
  CH = chow_ring(M)
  return [
    sum(
      t for t in Singular.terms(x) if __chow_ring_monomial_grading(M, t) == i; init=CH(0)
    )
    for
    i in 0:n
  ]
end

function truncate(M::QuiverModuliSpace, x, n)
  comps = homogeneous_components(M, x)
  # TODO this may be slow, try building new polynomial
  # also can I use the incorrect but faster Singular.degree?
  return sum(comps[i + 1] for i in 0:n)
end

"""
    _chern_characters_wedge(F::Bundle, k)

Compute the exterior powers of `F` up to degree `k`.
For internal use only.
"""
function _chern_characters_wedge(F::Bundle, k)
  k == 0 && return [1]
  x = chern_character(F)
  M = variety(F)
  CH = chow_ring(F)
  n = dimension(M)

  # init as CH(0) for type stability
  wedges = [CH(0) for _ in 1:(k + 1)]
  wedges[1], wedges[2] = CH(1), x

  for j in 2:k
    wedges[j + 1] =
      CH(1//j) * truncate(M,
        sum(
          (-CH(1))^(j - i + 1) * wedges[i + 1] * adams(F, j - i) for i in 0:(j - 1);
          init=CH(0),
        ),
        n)
    simplify!(wedges[j + 1])
  end
  return wedges
end

"""
    _chern_characters_symmetric(F::Bundle, k)

Compute the symmetric powers of `F` up to degree `k`.
For internal use only.
"""
function _chern_characters_symmetric(F::Bundle, k)
  k == 0 && return [1]
  x = chern_character(F)
  M = variety(F)
  n = dimension(M)
  CH = chow_ring(F)
  r = rank(F)

  wedges = _chern_characters_wedge(F, r)
  # init as CH(0) for type stability
  syms = [CH(0) for _ in 1:(k + 1)]
  syms[1], syms[2] = CH(1), x

  for j in 2:k
    syms[j + 1] = truncate(M,
      sum(
        (-CH(1))^(i + 1) * wedges[i + 1] * syms[j - i + 1] for i in 1:min(j, r);
        init=CH(0),
      ),
      n)
    simplify!(syms[j + 1])
  end
  return syms
end

"""
    adams(F::Bundle, k)

Compute the Adams operation ``\\Phi^k`` on the Chern character of `F`.
For internal use only.
"""
function adams(F::Bundle, k)
  M = variety(F)
  n = dimension(M)
  x = chern_character(F)

  return [k^i for i in 0:n]' * homogeneous_components(M, x)
end

function _chern_classes_from_character(F::Bundle)
  CH = chow_ring(F)
  M = variety(F)
  n = dimension(M)
  comps = homogeneous_components(M, chern_character(F))
  p = [(CH(-1))^i * CH(factorial(i)) * comps[i + 1] for i in 0:n]
  e = [CH(0) for _ in 1:(n + 1)]
  e[1] = CH(1)
  for i in 1:n
    e[i + 1] = CH(-1//i) * sum(p[j + 1] * e[i - j + 1] for j in 1:i)
    simplify!(e[i + 1])
  end
  return Dict(i => e[i + 1] for i in 0:n)
end

function _chern_character_from_classes(F::Bundle)
  CH = F.parent.ring
  M = variety(F)
  n = dimension(M)
  n == 0 && return CH(0)
  e = chern_classes(F)
  p = vcat([-e[1]], [CH(0) for _ in 1:(n - 1)])
  for i in 1:(n - 1)
    p[i + 1] = -CH(i + 1) * e[i + 1] - sum(e[j] * p[i - j + 1] for j in 1:i)
  end
  return simplify(sum(CH((-1)^i//factorial(i)) * p[i] for i in 1:n) + rank(F))
end

simplify(f::Singular.spoly{Singular.n_Q}) = div(f, f.parent(1))

function simplify!(f::Singular.spoly{Singular.n_Q})
  f = div(f, f.parent(1))
  return f
end

"""
    canonical_bundle(M::QuiverModuliSpace)

Return the canonical bundle on the quiver moduli space `M`.

If ``d`` is ``theta``-coprime and amply stable, the canonical bundle
is described in [Proposition 4.2, MR4352662](https://mathscinet.ams.org/mathscinet-getitem?mr=4352662).

# Example

On various projective spaces:

```jldoctest
julia> d = [1, 1]; Pn = map(i -> QuiverModuliSpace(kronecker_quiver(i + 1), d), 1:5);

julia> omega = map(canonical_bundle, Pn);

julia> map(chern_class, omega)
5-element Vector{Singular.spoly{Singular.n_Q}}:
 2*x11
 3*x11
 4*x11
 5*x11
 6*x11
```
"""
function canonical_bundle(M::QuiverModuliSpace)
  !(is_coprime(M) && is_amply_stable(M)) &&
    throw(
      NotImplementedError(
        "not coprime and amply stable, cannot compute the canonical bundle."
      ),
    )
  cl_omega = chern_class_line_bundle(M, -canonical_stability(M.Q, M.d))
  new = Bundle(M, 1, cl_omega)

  r = index(M)
  weights = all_weights_irreducible_component_canonical(M)
  for hn in keys(weights)
    weights[hn] *= r
  end
  return set_bundle_weights!(new, weights)
end

"""
    universal_bundle(M:::QuiverModuliSpace, i::Int)

Returns the `i`-th universal bundle of `M`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=[1, 0]);

julia> u1, u2 = universal_bundle(M, 1), universal_bundle(M, 2);

julia> map(chern_class, [u1, u2])
2-element Vector{Singular.spoly{Singular.n_Q}}:
 x11 + 1
 x21 + 1
```
"""
function universal_bundle(M::QuiverModuliSpace, i::Int)
  !is_coprime(M) &&
    throw(
      ArgumentError("$(M.d) is not $(M.theta)-coprime, universal bundles do not exist.")
    )
  cl = total_chern_class_universal(M, i)
  new = Bundle(M, M.d[i], cl)

  weights = all_weights_universal_bundle(M, i)
  return set_bundle_weights!(new, weights)
end

"""
    degree(F::Bundle)

Return the degree of the bundle `F`.
If `rank(F)` is larger than ``1``, returns the degree of the determinant of `F`.

# Example

The degrees of canonical bundles on the first projective spaces:

```jldoctest
julia> d = [1, 1]; Pn = map(i -> QuiverModuliSpace(kronecker_quiver(i + 1), d), 1:5);

julia> omega = map(canonical_bundle, Pn);

julia> omega = map(dual, omega);

julia> map(degree, omega)
5-element Vector{Singular.spoly{Singular.n_Q}}:
 2
 9
 64
 625
 7776
```

An example from [arXiv:2411.15125](https://arxiv.org/abs/2411.15125):

```jldoctest
julia> Q = Quiver("1-2,1-3,2---3"); d = [1, 1, 1]; a = [1, 1, -1];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> F = dual(canonical_bundle(M))
Bundle of rank 1.

julia> chern_class(F)
2*x31 + 1

julia> degree(F)
56
"""
function degree(F::Bundle)
  M = variety(F)
  n = dimension(M)
  return div(homogeneous_components(M, chern_class(det(F))^n)[n + 1], point_class(M))
end
