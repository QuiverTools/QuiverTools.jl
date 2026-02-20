###############################################
# Methods related to stability and lack thereof
###############################################

"""
    canonical_stability(Q::Quiver, d::AbstractVector{Int})

Compute the canonical stability parameter for `Q` and `d`.

This is defined to be ``<d,-> - <-,d>``

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.

# Output

- the canonical stability parameter for `Q` and `d`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2,3];

julia> canonical_stability(Q, d) == [9, -6]
true
```
"""
function canonical_stability(Q::Quiver, d::AbstractVector{Int})
  return -(-transpose(euler_matrix(Q)) + euler_matrix(Q)) * d
end

"""
    is_coprime(d::AbstractVector{Int}, theta::AbstractVector{Int})

Check if `d` is `theta`-coprime.

A dimension vector ``d`` is said to be ``\\theta``-coprime for the
stability parameter ``\\theta`` if all subdimension vectors ``0 \\neq e < d``
satisfy ``\\theta * e \\neq 0``.

# Input

- `d::AbstractVector{Int}` a dimension vector.
- `theta::AbstractVector{Int}` a stability parameter.

# Output

- `true` if `d` is `theta`-coprime, `false` otherwise.

# Examples

```jldoctest
julia> d = [2, 3]; theta = [3, -2];

julia> is_coprime(d, theta)
true

julia> is_coprime([3, 3], theta)
false
```
"""
function is_coprime(d::AbstractVector{Int}, theta::AbstractVector{Int})
  return all(
    e -> theta' * e != 0,
    all_subdimension_vectors(d; nonzero=true, strict=true),
  )
end

"""
    is_coprime(d::AbstractVector{Int})

Check if the gcd of all the entries of d is ``1``.

# Input

- `d::AbstractVector{Int}` a vector.

# Output

- `true` if the gcd of all the entries of `d` is ``1``, `false` otherwise.

# Examples

```jldoctest
julia> is_coprime([2, 3])
true

julia> is_coprime([3, 3])
false
```
"""
is_coprime(d::AbstractVector{Int}) = gcd(d) == 1

"""
    slope(d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)

Return the slope of `d`
with respect to the stability parameter `theta`
and a choice of a denominator function `denom`.

The slope function for ``\\theta`` and ``\\alpha`` is defined as
```math
\\mu = \\frac{\\theta}{\\alpha} := x \\mapsto \\frac{\\theta \\cdot x}{\\alpha(x)}.
```

# Input

- `d::AbstractVector{Int}` a dimension vector.
- `theta::AbstractVector{Int}` a stability parameter.
- `denom::Function` a function to compute the denominator. Default is `sum`.

# Output

- the slope of `d` with respect to `theta` and `denom`.

# Examples

```jldoctest
julia> slope([2,3], [3,-2])
0//1
```
"""
function slope(d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
  return (theta' * d)//denom(d)
end

"""
    all_destabilizing_subdimension_vectors(d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)

Return the subdimension vectors of `d` with a strictly larger slope than `d`.

# Input

- `d::AbstractVector{Int}` a dimension vector.
- `theta::AbstractVector{Int}` a stability parameter.
- `denom::Function` a function to compute the denominator. Default is `sum`.

# Output

- an array of subdimension vectors of `d` with a strictly larger slope than `d`.

# Examples

```jldoctest
julia> QuiverTools.all_destabilizing_subdimension_vectors([2, 3], [3, -2])
5-element Vector{Vector{Int64}}:
 [1, 0]
 [2, 0]
 [1, 1]
 [2, 1]
 [2, 2]

julia> QuiverTools.all_destabilizing_subdimension_vectors([2, 3], [0, 0])
Vector{Int64}[]

julia> QuiverTools.all_destabilizing_subdimension_vectors([0, 0], [1, -1])
Vector{Int64}[]
```
"""
function all_destabilizing_subdimension_vectors(
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  all(di == 0 for di in d) && return Vector{Int}[]

  b = slope(d, theta, denom)
  return filter(
    e -> slope(e, theta, denom) > b,
    all_subdimension_vectors(d; nonzero=true, strict=true),
  )
end

"""
    has_semistables(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d), denom::Function=sum)

Check if `Q` admits a `theta`-semistable representation of dimension vector `d`.

A representation ``V`` is said to be ``\\mu``-semistable if
for all of its subrepresentations ``W`` with ``dim(W) < dim(V)``, we have
```math
\\mu(\\dim(W)) \\leq \\mu\\dim((V)).
```

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.
- `theta::AbstractVector{Int}` a stability parameter. Default is `canonical_stability(Q, d)`.
- `denom::Function` a function to compute the denominator. Default is `sum`.

# Output

- `true` if there is a `theta`-semistable representation of dimension vector `d`,
  `false` otherwise.

# Examples

```jldoctest
julia> A2 = kronecker_quiver(1); theta = [1,-1];

julia> has_semistables(A2, [1,1], theta)
true

julia> has_semistables(A2, [2,2], theta)
true

julia> has_semistables(A2, [1,2], theta)
false

julia> has_semistables(A2, [0,0], theta)
true
```
The 3-Kronecker quiver:

```jldoctest
julia> K3 = kronecker_quiver(3); theta = [3,-2];

julia> has_semistables(K3, [2,3], theta)
true

julia> has_semistables(K3, [1,4], theta)
false
```
"""
@memoize Dict function has_semistables(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
)
  all(di == 0 for di in d) && return true

  # collect the list of all subdimension vectors e of bigger slope than d
  slope_d = slope(d, theta, denom)

  subdimensions_bigger_slope = filter(
    e -> slope(e, theta, denom) > slope_d,
    all_subdimension_vectors(d; nonzero=true, strict=true),
  )
  # to have semistable representations, none of the vectors above must be
  # a general subdimension vector.
  return all(e -> !is_general_subdimension_vector(Q, e, d), subdimensions_bigger_slope)
end

"""
    has_stables(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d), denom::Function=sum)

Check if Q admits a `theta`-stable representation of dimension vector `d`.

A representation ``V`` is said to be ``mu``-stable if
for all of its subrepresentations ``W`` with ``dim(W) < dim(V)``, we have
```math
\\mu(\\dim(W)) < \\mu\\dim((V)).
```

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.
- `theta::AbstractVector{Int}` a stability parameter. Default is `canonical_stability(Q, d)`.
- `denom::Function` a function to compute the denominator. Default is `sum`.

# Output

- `true` if there is a `theta`-stable representation of dimension vector `d`,
  `false` otherwise.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2, 3]; theta = [3, -2];

julia> has_stables(Q, d, theta)
true

julia> Q = kronecker_quiver(2); d = [2,2]; theta = [1,-1];

julia> has_stables(Q, d, theta)
false

julia> has_semistables(Q, d, theta)
true
```

The zero dimension vector has no stables:
```jldoctest
julia> Q = kronecker_quiver(3); d = [0,0];

julia> has_stables(Q, d)
false
```
"""
@memoize Dict function has_stables(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
)
  all(di == 0 for di in d) && return false
  # collect the list of all subdimension vectors e of bigger slope than d
  slope_d = slope(d, theta, denom)

  subdimensions_bigger_or_equal_slope = filter(
    e -> slope(e, theta, denom) >= slope_d,
    all_subdimension_vectors(d; nonzero=true, strict=true),
  )
  # to have semistable representations,
  # none of the vectors above must be general subdimension vectors.
  return all(
    e -> !is_general_subdimension_vector(Q, e, d),
    subdimensions_bigger_or_equal_slope,
  )
end

# TODO the cited paper is published
# The published version is not open access.
"""
    is_schur_root(Q::Quiver, d::AbstractVector{Int})

Check if `d` is a Schur root for `Q`.

By [[Lemma 4.2, arXiv:0802.2147](https://doi.org/10.48550/arXiv.0802.2147)],
this is equivalent to the existence of a stable representation of dimension vector ``d``
for the canonical stability parameter.

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.

# Output

- `true` if `d` is a Schur root for `Q`, `false` otherwise.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2,3];

julia> is_schur_root(Q, d)
true
```
"""
is_schur_root(Q::Quiver, d::AbstractVector{Int}) = has_stables(
  Q, d, canonical_stability(Q, d)
)

"""
    is_root(Q::Quiver, d)

Check whether `d` is a root, i.e., if ``<d, d> \\leq 1``.
"""
is_root(Q::Quiver, d) = euler_form(Q, d, d) <= 1

"""
    is_real_root(Q::Quiver, d)

Check whether `d` is a real root, i.e., if ``<d, d> = 1``.
"""
is_real_root(Q::Quiver, d) = euler_form(Q, d, d) == 1

"""
    is_imaginary_root(Q::Quiver, d)

Check whether `d` is an imaginary root, i.e., if ``<d, d> \\geq 0``.
"""
is_imaginary_root(Q::Quiver, d) = euler_form(Q, d, d) <= 0

"""
    is_isotropic_root(Q::Quiver, d)

Check whether `d` is an isotropic root, i.e., if ``<d, d> = 0``.
"""
is_isotropic_root(Q::Quiver, d) = euler_form(Q, d, d) == 0

"""
    is_general_subdimension_vector(Q::Quiver, e::AbstractVector{Int}, d::AbstractVector{Int})

Check if `e` is a general subdimension vector of `d`.

A dimension vector ``e`` is called a general subdimension vector of ``d``
if a general representation of dimension vector ``d`` possesses a subrepresentation
of dimension vector ``e``.

By [[Theorem 5.3, arXiv:0802.2147](https://doi.org/10.48550/arXiv.0802.2147)],
``e`` is a general subdimension vector of ``d`` if and only if
```math
<e',d-e> \\geq 0
```
for all general subdimension vectors ``e'`` of ``e``.

# Input

- `Q::Quiver` a quiver.
- `e::AbstractVector{Int}` a dimension vector.
- `d::AbstractVector{Int}` a dimension vector.

# Output

- `true` if `e` is a general subdimension vector of `d`, `false` otherwise.

# Examples

Trivial examples on the 3-Kronecker quiver:

```jldoctest
julia> Q = kronecker_quiver(3); e = [1, 2]; d = [2, 3];

julia> is_general_subdimension_vector(Q, e, d)
true

julia> is_general_subdimension_vector(Q, [0, 0], d)
true

julia> is_general_subdimension_vector(Q, [2, 3], d)
true

julia> is_general_subdimension_vector(Q, [2, 1], d)
false
```
"""
@memoize Dict function is_general_subdimension_vector(
  Q::Quiver,
  e::AbstractVector{Int},
  d::AbstractVector{Int},
)
  (e == d || all(ei == 0 for ei in e)) && return true

  # to speed up computation of <eprime,d-e>
  partial_evaluation = euler_matrix(Q) * (d - e)
  # considering subdimension vectors that violate the numerical condition
  subdimensions = filter(
    eprime -> eprime' * partial_evaluation < 0, all_subdimension_vectors(e)
  )
  # none of the subdimension vectors violating the condition should be general
  return all(eprime -> !is_general_subdimension_vector(Q, eprime, e), subdimensions)
  # return general_ext(Q, e, d - e) == 0 # TODO test performance
end

"""
    all_general_subdimension_vectors(Q::Quiver, d::AbstractVector{Int})

Return the list of all general subdimension vectors of `d`.

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.

# Output

- a list of all general subdimension vectors of `d`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> QuiverTools.all_general_subdimension_vectors(Q, [2, 3])
7-element Vector{Vector{Int64}}:
 [0, 0]
 [0, 1]
 [0, 2]
 [1, 2]
 [0, 3]
 [1, 3]
 [2, 3]

julia> QuiverTools.all_general_subdimension_vectors(Q, [3, 0])
4-element Vector{Vector{Int64}}:
 [0, 0]
 [1, 0]
 [2, 0]
 [3, 0]
```
"""
@memoize Dict function all_general_subdimension_vectors(Q::Quiver, d::AbstractVector{Int})
  return filter(e -> is_general_subdimension_vector(Q, e, d), all_subdimension_vectors(d))
end

"""
    all_hn_types(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum; unstable::Bool=false, ordered::Bool=true)

Return a list of all the Harder--Narasimhan types of representations of `Q`
with dimension vector `d`, with respect to the slope function `theta`/`denom`.

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.
- `theta::AbstractVector{Int}` a stability parameter.
- `denom::Function` a function to compute the denominator. Default is `sum`.

Keyword inputs:

- `unstable`: if `true` exclude the trivial Harder--Narasimhan type (d),
which corresponds to stable representations. Default is `false`.
- `ordered`: if `true` return the list of all Harder--Narasimhan types in ascending order.
Default is `true`.

# Output

- a list of all the Harder--Narasimhan types
of representations of `Q` with dimension vector `d`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> all_hn_types(Q, d, theta; ordered=true)
8-element Vector{HNType}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> all_hn_types(Q, [3,0], [0,0]) == [[[3, 0]]]
true

julia> Q = three_vertex_quiver(1, 4, 1); d = [4, 1, 4];

julia> theta = canonical_stability(Q, d);

julia> length(all_hn_types(Q, d, theta; ordered=true))
106
```
"""
@memoize Dict function all_hn_types(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum;
  unstable::Bool=false,
  ordered::Bool=true,
)
  all(di == 0 for di in d) && return [HNType([zero_vector(Q)])]

  # We consider just proper subdimension vectors which admit a semistable
  # representation and for which μ(e) > μ(d)
  # Note that we also eliminate d by the following
  subdimensions = all_destabilizing_subdimension_vectors(d, theta, denom)
  filter!(e -> has_semistables(Q, e, theta, denom), subdimensions)

  # We sort the subdimension vectors by slope because that will return the list of
  # all HN types in ascending order with respect to the partial order from
  # Definition 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
  ordered && sort!(subdimensions; by=e -> slope(e, theta, denom))

  # The HN types which are not of the form (d) are (e,f^1,...,f^s) where e is a
  # proper semistable subdimension vector with μ(e) > μ(d), (f^1,...,f^s) is a HN
  # type of f = d-e and μ(e) > μ(f^1) holds.

  alltypes = HNType[
    HNType(vcat([e], efstar.hn))

    for e in subdimensions for efstar in filter(
      fstar -> slope(e, theta, denom) > slope(fstar[1], theta, denom),
      all_hn_types(Q, d - e, theta, denom; ordered=ordered),
    )
  ]

  # Possibly add d again, at the beginning, because it is smallest
  # with respect to the partial order from Definition 3.6
  if !unstable && has_semistables(Q, d, theta, denom)
    pushfirst!(alltypes, HNType([d]))
  end
  return alltypes
end

"""
    is_hn_type(Q::Quiver, d::AbstractVector{Int}, dstar::HNType, theta::AbstractVector{Int}=canonical_stability(Q, d), denom::Function=sum)

Check if the given ordered list of subdimension vectors `dstar` is a Harder--Narasimhan type
for the datum `Q`, `d` and the slope function `theta`/`denom`.

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.
- `dstar::HNType` an Harder--Narasimhan type.
- `theta::AbstractVector{Int}` a stability parameter.
- `denom::Function` a function to compute the denominator. Default is `sum`.

# Output

- `true` if `dstar` is a Harder--Narasimhan type for `Q`, `d` and the slope `theta`/`denom`, `false` otherwise.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> d = [2, 3]; dstar = [d];

julia> is_hn_type(Q, d, dstar)
true
```
"""
function is_hn_type(
  Q::Quiver,
  d::AbstractVector{Int},
  dstar::HNType,
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
)
  sum(dstar) != d && throw(ArgumentError("$(dstar) does not sum to $(d)."))

  if !all(
    slope(dstar[i], theta, denom) > slope(dstar[i + 1], theta, denom) for
    i in 1:(length(dstar) - 1)
  )
    return false
  end

  if !all(has_semistables(Q, dstari, theta, denom) for dstari in dstar)
    return false
  end
  return true
end
is_hn_type(Q::Quiver,
d::AbstractVector{Int},
dstar::Vector{<:AbstractVector{Int}};
theta::AbstractVector{Int}=canonical_stability(Q, d),
denom::Function=sum
) = is_hn_type(Q, d, HNType(dstar), theta, denom)

"""
    codimension_hn_stratum(Q::Quiver, stratum::HNType)

Compute the codimension of the given Harder--Narasimhan stratum.

# Input

- `Q::Quiver` a quiver.
- `stratum::HNType`: a Harder--Narasimhan type.

# Output

- the codimension of the Harder--Narasimhan stratum as an integer.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> HN = all_hn_types(Q, d, theta; ordered=true);

julia> [codimension_hn_stratum(Q, stratum) for stratum in HN]
8-element Vector{Int64}:
  0
  3
  4
 10
  8
  9
 12
 18
```
"""
function codimension_hn_stratum(Q::Quiver, stratum::HNType)
  length(stratum) == 1 && return 0

  return -sum(
    euler_form(Q, stratum[i], stratum[j])
    for i in 1:(length(stratum) - 1)
    for j in (i + 1):length(stratum); init=0
  )
end
codimension_hn_stratum(Q::Quiver, stratum::Vector{<:AbstractVector{Int}}) = codimension_hn_stratum(
  Q, HNType(stratum)
)

"""
    is_amply_stable(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)

Check whether the dimension vector `d` is amply stable
with respect to the slope function `theta`/`denominator`.

This means that the codimension of the unstable locus
in the parameter space is at least ``2``.

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.
- `theta::AbstractVector{Int}` a stability parameter.
- `denom::Function` a function to compute the denominator. Default is `sum`.

# Output

- `true` if `d` is amply stable, `false` otherwise.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2, 3];

julia> is_amply_stable(Q, d, [3, -2])
true

julia> is_amply_stable(Q, d, [-3, 2])
false

julia> is_amply_stable(Q, [3, 0], [0, -3])
true
```
"""
function is_amply_stable(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
)
  hn_types = all_hn_types(Q, d, theta; unstable=true, ordered=false)
  return all(stratum -> codimension_hn_stratum(Q, stratum) >= 2, hn_types)
end

"""
    has_properly_semistables(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, denom::Function=sum)

"""
function has_properly_semistables(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum
)
  is_coprime(d, theta) && return false
  return !isempty(all_luna_types(Q, d, theta, denom; stable=false))
end
