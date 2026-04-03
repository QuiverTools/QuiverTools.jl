######################################################################
# Weights of various standard vector bundles for the HN stratification
######################################################################

"""
    weights_hn_type(hntype::HNType, theta::AbstractVector{Int}, denom::Function=sum)

Compute the weights of the 1-PS corresponding to `hn_type` for the slope function
`theta`/`denom`.

# Example

For our favourite 6-fold, the weights of the 1-PS corresponding to each HN type are:

```jldoctest
julia> Q = kronecker_quiver(3); d = [2, 3]; theta = [3, -2];

julia> hn = all_hn_types(Q, d, theta; unstable=true);

julia> map(hn_type -> weights_hn_type(hn_type, theta), hn)
7-element Vector{Vector{Int64}}:
 [3, -2]
 [1, -4]
 [2, -3]
 [4, -1]
 [9, -1, -6]
 [6, 1, -4]
 [3, -2]
```
"""
function weights_hn_type(hntype::HNType, theta::AbstractVector{Int}, denom::Function=sum)
  k_weights = map(h -> slope(h, theta, denom), hntype)

  c = lcm(denominator.(k_weights))
  map!(k -> c * k, k_weights, k_weights)
  gg = gcd(k_weights)
  map!(k -> k / gg, k_weights, k_weights)
  return map(Int, k_weights) # if we accepted rational weights this would not be needed.
end

"""
    weight_line_bundle_on_stratum(hn_type::HNType, eta::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)

Compute the weight on `hn_type` of the line bundle with linearization `eta`.
"""
function weight_line_bundle_on_stratum(
  hn_type::HNType,
  eta::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  k_weights = weights_hn_type(hn_type, theta, denom)
  return [
    Int(
      -eta' * sum(
        k_weights[m] .* hn_type[m] for m in 1:length(hn_type);
        init=zeros(Int, length(hn_type[1])),
      ),
    ),
  ]
end

"""
    weights_line_bundle(Q::Quiver, d::AbstractVector{Int}, eta::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)

Compute the Teleman weights on the line bundle given by the linearization `eta`.
"""
function weights_line_bundle(Q::Quiver,
  d::AbstractVector{Int},
  eta::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  hn_types = all_hn_types(Q, d, theta, denom; unstable=true)
  return Dict(
    hn_type => weight_line_bundle_on_stratum(hn_type, eta, theta, denom) for
    hn_type in hn_types
  )
end

function weights_line_bundle(M::QuiverModuliSpace, eta::AbstractVector{Int})
  hn_types = all_hn_types(M.Q, M.d, M.theta, M.denom; unstable=true)
  return Dict(
    hn_type => weight_line_bundle_on_stratum(hn_type, eta, M.theta, M.denom) for
    hn_type in hn_types
  )
end

"""
    teleman_bound_on_stratum(Q::Quiver, hn_type::HNType, theta::AbstractVector{Int}, denom::Function=sum)

Compute the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS ``\\lambda``
corresponding to the given HN type.

# Input

- `Q`: a quiver.
- `hn_type`: a Harder--Narasimhan type.
- `theta`: a stability parameter.
- `denom`: a denominator for the slope function. Defaults to `sum`.

# Output

The weight of the 1-PS corresponding to the given HN type.
"""
function teleman_bound_on_stratum(
  Q::Quiver,
  hn_type::HNType,
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  ell = length(hn_type)
  ell == 1 &&
    throw(ArgumentError("Weight not defined on the dense stratum"))

  k_weights = weights_hn_type(hn_type, theta, denom)
  return sum(
    (k_weights[t] - k_weights[s]) * euler_form(Q, hn_type[s], hn_type[t])
    for s in 1:(ell - 1) for t in (s + 1):ell; init=0
  )
end

function teleman_bound_on_stratum(M::QuiverModuli, hn_type::HNType)
  return teleman_bound_on_stratum(M.Q, hn_type, M.theta, M.denom)
end

"""
    teleman_bounds(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)

Compute the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each
HN type for the given `Q`, `d`, `\\theta` and `denom`.

# Input

- `Q::Quiver`: a quiver
- `d::AbstractVector{Int}`: a dimension vector
- `theta::AbstractVector{Int}`: a stability parameter
- `denom::Function`: a denominator for the slope function. Defaults to `sum`.

# Output

A dictionary with the weights of the 1-PS corresponding to each HN type.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> teleman_bounds(Q, [2, 3], [3, -2])
Dict{HNType, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 20
  [[2, 1], [0, 2]]         => 50
  [[1, 0], [1, 2], [0, 1]] => 100
  [[1, 0], [1, 3]]         => 40
  [[1, 0], [1, 1], [0, 2]] => 90
  [[1, 1], [1, 2]]         => 15
  [[2, 0], [0, 3]]         => 90
```
"""
function teleman_bounds(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  hn_types = all_hn_types(Q, d, theta, denom; unstable=true)
  return Dict(
    hn_type => teleman_bound_on_stratum(Q, hn_type, theta, denom) for hn_type in hn_types
  )
end

"""
    teleman_bounds(M::QuiverModuli)

Compute the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each
HN type for the datum of `M`.

# Input

- `M::QuiverModuli`: a quiver moduli space or stack

# Output

A dictionary with the weights of the 1-PS corresponding to each HN type.

# Examples

```jldoctest
julia> Q = three_vertex_quiver(1, 2, 3); d = [3, 1, 2]; theta = [5, 3, -9];

julia> M = QuiverModuliSpace(Q, d, theta);

julia> length(teleman_bounds(M))
24
```
"""
function teleman_bounds(M::QuiverModuli)
  return teleman_bounds(M.Q, M.d, M.theta, M.denom)
end

"""
    weights_universal_bundle_on_stratum(hn_type::HNType, i::Int, theta::AbstractVector{Int}, denom::Function=sum; chi::AbstractVector{Int})

Returns the weights of a universal bundle ``U_i(a)`` for the linearization ``a``
for the 1-PS corresponding to the given HN type.

"""
function weights_universal_bundle_on_stratum(
  hn_type::HNType,
  i::Int,
  theta::AbstractVector{Int},
  denom::Function=sum;
  chi::AbstractVector{Int},
)
  ell = length(hn_type)
  k_weights = weights_hn_type(hn_type, theta, denom)
  constant_term = sum(k_weights[s] * (chi' * hn_type[s]) for s in 1:ell; init=0)

  weights_mult = reduce(
    vcat, [k_weights[s] for _ in 1:hn_type[s][i]] for s in 1:ell
  )
  map!(w -> -constant_term + w, weights_mult, weights_mult)
  return weights_mult
end

"""
    weights_universal_bundle(Q::Quiver, d::AbstractVector{Int}, i::Int, theta::AbstractVector{Int}, denom::Function=sum; chi::AbstractVector{Int})

Compute the Teleman weights of the universal bundle ``U_i(chi)``
for the linearization ``chi`` on all the non-dense Harder-Narasimhan strata.

# Input

- `Q::Quiver`: a quiver
- `d::AbstractVector{Int}`: a dimension vector
- `i::Int`: the index of the universal bundle
- `theta::AbstractVector{Int}`: a stability parameter
- `denom::Function`: a denominator for the slope function. Defaults to `sum`.

Keyword arguments:

- `chi::AbstractVector{Int}`: the linearization of the universal bundle. Defaults to `extended_gcd(d)[2]`.

# Output

A dictionary with the weights of the universal bundle on each stratum.
"""
function weights_universal_bundle(
  Q::Quiver,
  d::AbstractVector{Int},
  i::Int,
  theta::AbstractVector{Int},
  denom::Function=sum;
  chi::AbstractVector{Int},
)
  gcd(d) > 1 &&
    throw(
      ArgumentError(
        "gcd($(collect(d))) = $(gcd(d)) > 1, the universal bundles do not exist."
      ),
    )
  hn_types = all_hn_types(Q, d, theta, denom; unstable=true)
  return Dict(
    hn_type => weights_universal_bundle_on_stratum(hn_type, i, theta, denom; chi=chi)
    for hn_type in hn_types
  )
end

"""
    weights_universal_bundle(M::QuiverModuli, i::Int; chi::Union{AbstractVector{Int},UndefInitializer}=undef)

Compute the Teleman weights of the universal bundle ``U_i(chi)``
for the linearization `chi` on all the non-dense Harder-Narasimhan strata.

# Input

- `M::QuiverModuli`: a quiver moduli space or stack.

Keyword arguments:

- `chi::AbstractVector{Int}`: the linearization of the universal bundle. Defaults to `linearization(M)`.

# Output

A dictionary with the weights of the universal bundle ``U_{i}(chi)`` on each stratum.

# Example

The weights of the universal bundles on our favourite 6-fold:

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> weights_universal_bundle(M, 1; chi=[2, -1])
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [-5, -5]
  [[2, 1], [0, 2]]         => [-10, -10]
  [[1, 0], [1, 2], [0, 1]] => [-15, -25]
  [[1, 0], [1, 3]]         => [-5, -10]
  [[1, 0], [1, 1], [0, 2]] => [-15, -20]
  [[1, 1], [1, 2]]         => [0, -5]
  [[2, 0], [0, 3]]         => [-15, -15]
```

If not specified, the linearization is taken from the moduli space, and
defaults to `extended_gcd(M.d)[2]` if not defined.

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> weights_universal_bundle(M, 1)
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [5, 5]
  [[2, 1], [0, 2]]         => [10, 10]
  [[1, 0], [1, 2], [0, 1]] => [25, 15]
  [[1, 0], [1, 3]]         => [10, 5]
  [[1, 0], [1, 1], [0, 2]] => [20, 15]
  [[1, 1], [1, 2]]         => [5, 0]
  [[2, 0], [0, 3]]         => [15, 15]

julia> QuiverTools.set_linearization!(M, [-4, 3]);

julia> weights_universal_bundle(M, 1)
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [15, 15]
  [[2, 1], [0, 2]]         => [30, 30]
  [[1, 0], [1, 2], [0, 1]] => [65, 55]
  [[1, 0], [1, 3]]         => [25, 20]
  [[1, 0], [1, 1], [0, 2]] => [55, 50]
  [[1, 1], [1, 2]]         => [10, 5]
  [[2, 0], [0, 3]]         => [45, 45]
```
"""
function weights_universal_bundle(
  M::QuiverModuli,
  i::Int;
  chi::Union{AbstractVector{Int},UndefInitializer}=undef,
)
  if !(chi isa UndefInitializer)
    return weights_universal_bundle(M.Q, M.d, i, M.theta, M.denom; chi=chi)
  end

  if M isa QuiverModuliSpace
    return weights_universal_bundle(M.Q, M.d, i, M.theta, M.denom; chi=linearization(M))
  end

  gcd_value, default_chi = extended_gcd(M.d)
  gcd_value != 1 && throw(
    ArgumentError(
      "No default linearization exists because gcd($(collect(M.d))) = $(gcd_value) != 1."
    ),
  )
  return weights_universal_bundle(M.Q, M.d, i, M.theta, M.denom; chi=default_chi)
end

# TODO implement irreducible component as well.
"""
    weight_canonical_on_stratum(Q::Quiver, d::AbstractVector{Int}, hn_type::HNType, theta::AbstractVector{Int}, denom::Function=sum)

Compute the Teleman weight of ``\\omega_R|_Z``
on the Harder-Narasimhan stratum `hn_type`.

"""
function weight_canonical_on_stratum(
  Q::Quiver,
  d::AbstractVector{Int},
  hn_type::HNType,
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  return weight_line_bundle_on_stratum(hn_type, -canonical_stability(Q, d), theta, denom)
end

"""
    weight_irreducible_component_canonical_on_stratum(M::QuiverModuli, hn_type::HNType)

Compute the Teleman weight of the irreducible component of ``\\omega_R|_Z``
on the Harder-Narasimhan stratum `hn_type`.

More explicitly, if ``\\omega_X = \\mathcal{O}(rH)``, this returns the weight of
the pullback of O(H) on the given stratum.

"""
function weight_canonical_on_stratum(
  M::QuiverModuli,
  hn_type::HNType,
)
  return weight_canonical_on_stratum(
    M.Q,
    M.d,
    hn_type,
    M.theta,
    M.denom,
  )
end

"""
    weights_canonical_bundle(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)

Compute the Teleman weights of the irreducible component of ``\\omega_R|_Z``
on all the non-dense Harder-Narasimhan strata.

More explicitly, if ``\\omega_X = O(rH)``, this returns the weights of the pullback of
``\\mathcal{O}(H)`` on each stratum.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter.
- `denom::Function`: a denominator for the slope function. Defaults to `sum`.

# Output

A dictionary with the weights of the irreducible component of the canonical bundle on each stratum.
"""
function weights_canonical_bundle(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  has_properly_semistables(Q, d, theta, denom) &&
    throw(
      ArgumentError(
        "The quiver moduli problem has properly semistables, no description of the canonical bundle is available."
      ),
    )
  !is_amply_stable(Q, d, theta) &&
    throw(
      ArgumentError(
        "The quiver moduli problem is not amply stable, no description of the canonical bundle is available."
      ),
    )
  return weights_line_bundle(Q, d, -canonical_stability(Q, d), theta, denom)
end

"""
    weights_canonical_bundle(M::QuiverModuliSpace)

Compute the Teleman weights of ``\\omega_R|_Z``
on all the non-dense Harder-Narasimhan strata.

# Input

- `M::QuiverModuliSpace`: a quiver moduli space.

# Output

A dictionary with the weights of ``\\omega_{R}|_{Z}`` on each stratum ``Z``.

# Example

The canonical bundle of our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> weights_canonical_bundle(M)
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [30]
  [[2, 1], [0, 2]]         => [60]
  [[1, 0], [1, 2], [0, 1]] => [120]
  [[1, 0], [1, 3]]         => [45]
  [[1, 0], [1, 1], [0, 2]] => [105]
  [[1, 1], [1, 2]]         => [15]
  [[2, 0], [0, 3]]         => [90]
```
"""
function weights_canonical_bundle(M::QuiverModuliSpace)
  return weights_canonical_bundle(M.Q, M.d, M.theta, M.denom)
end

########################
# methods below shall be obsolete

"""
    weights_endomorphism_universal_bundle_on_stratum(hn_type::HNType, theta::AbstractVector{Int}, denom::Function=sum)

Compute all the weights that can occur in ``U^{\\vee} \\otimes U``
on the given Harder-Narasimhan stratum.

"""
function weights_endomorphism_universal_bundle_on_stratum(
  hn_type::HNType,
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  k_weights = weights_hn_type(hn_type, theta, denom)
  return [
    Int(k_weights[i] - k_weights[j]) for i in 1:length(hn_type) for j in 1:length(hn_type)
  ]
end

"""
    weights_endomorphism_universal_bundle_on_stratum(M::QuiverModuli, hn_type::HNType)

Compute all the weights that can occur in ``U^{\\vee} \\otimes U``
on the given Harder-Narasimhan stratum.

"""
function weights_endomorphism_universal_bundle_on_stratum(
  M::QuiverModuli,
  hn_type::HNType,
)
  return weights_endomorphism_universal_bundle_on_stratum(hn_type, M.theta, M.denom)
end

"""
    all_weights_endomorphisms_universal_bundle(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)

Compute all the weights that can occur in ``U^{\\vee} \\otimes U``
on the given Harder-Narasimhan stratum.

"""
function all_weights_endomorphisms_universal_bundle(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  gcd(d) > 1 && throw(
    ArgumentError("gcd($(M.d))  = $(gcd(d)) > 1, the universal bundles do not exist.")
  )

  hn_types = all_hn_types(Q, d, theta, denom; unstable=true)
  return Dict(
    hn_type => weights_endomorphism_universal_bundle_on_stratum(hn_type, theta, denom) for
    hn_type in hn_types
  )
end

"""
    all_weights_endomorphisms_universal_bundle(M::QuiverModuli)

Compute all the possible weights of the endomorphisms of the universal bundles
``U_i \\otimes U_j`` on all the non-dense Harder-Narasimhan strata
for each 1-PS relative to the HN type.

# Example

The weights of the endomorphisms of the universal bundles on our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> all_weights_endomorphisms_universal_bundle(M)
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [0, 5, -5, 0]
  [[2, 1], [0, 2]]         => [0, 5, -5, 0]
  [[1, 0], [1, 2], [0, 1]] => [0, 10, 15, -10, 0, 5, -15, -5, 0]
  [[1, 0], [1, 3]]         => [0, 5, -5, 0]
  [[1, 0], [1, 1], [0, 2]] => [0, 5, 10, -5, 0, 5, -10, -5, 0]
  [[1, 1], [1, 2]]         => [0, 5, -5, 0]
  [[2, 0], [0, 3]]         => [0, 5, -5, 0]
```
"""
function all_weights_endomorphisms_universal_bundle(M::QuiverModuli)
  return all_weights_endomorphisms_universal_bundle(M.Q, M.d, M.theta, M.denom)
end

# methods above shall be obsolete
##########################

#####################################################
# These return the correct multiplicities for weights.
# They can be used when Teleman fails for all of the pairs, but still holds for most.
#####################################################

"""
    weights_endomorphisms_universal_bundles_on_stratum(i::Int, j::Int, hn_type::HNType, theta::AbstractVector{Int}, denom::Function=sum)

Compute the weights of ``U_i^{\\vee} \\otimes U_j`` on `hn_type` for the slope `theta`/`denom`,
with the correct multiplicities.
"""
function weights_endomorphisms_universal_bundles_on_stratum(
  i::Int,
  j::Int,
  hn_type::HNType,
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  ell = length(hn_type)
  k_weights = weights_hn_type(hn_type, theta, denom)

  return reduce(
    vcat,
    Int[k_weights[t] - k_weights[s] for _ in 1:(hn_type[s][i] * hn_type[t][j])]
    for s in 1:ell for t in 1:ell
  )
end

"""
    weights_endomorphisms_universal_bundles(Q::Quiver, d::AbstractVector{Int}, i::Int, j::Int, theta::AbstractVector{Int}, denom::Function=sum)
"""
function weights_endomorphisms_universal_bundles(
  Q::Quiver,
  d::AbstractVector{Int},
  i::Int,
  j::Int,
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  hn_types = all_hn_types(Q, d, theta, denom; unstable=true)
  return Dict(
    hn_type =>
      weights_endomorphisms_universal_bundles_on_stratum(i, j, hn_type, theta, denom)
    for hn_type in hn_types
  )
end

"""
    weights_endomorphisms_universal_bundles(M::QuiverModuli, i::Int, j::Int)
"""
function weights_endomorphisms_universal_bundles(
  M::QuiverModuli,
  i::Int,
  j::Int,
)
  return weights_endomorphisms_universal_bundles(
    M.Q, M.d, i, j, M.theta, M.denom
  )
end

"""
    does_rigidity_inequality_hold(M::QuiverModuli)

Check if the Teleman quantization criterion of
[[arXiv:2311.17003](https://doi.org/10.48550/arXiv.2311.17003)] holds.

In case the quiver is acyclic, this ensures that the moduli space is infinitesimally
rigid.

# Input

- `M::QuiverModuli`: a quiver moduli space or stack.

# Output

- `true` if the Teleman inequality holds, `false` otherwise.

# Examples

Our favourite 6-fold is rigid:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> does_rigidity_inequality_hold(M)
true
```

Some quiver moduli are rigid, but it can't be proved by this criterion:
the following moduli space can be shown to be ``\\mathbb{P}^6``, whose rigidity follows
from the Euler sequence (see
[[Example 4.8, arXiv:2311.17003](https://doi.org/10.48550/arXiv.2311.17003)]). However,
the Teleman inequality does not hold:
```jldoctest
julia> Q = three_vertex_quiver(1, 6, 1); M = QuiverModuliSpace(Q, [1, 6, 6], [42, 5, -12]);

julia> does_rigidity_inequality_hold(M)
false
```
"""
function does_rigidity_inequality_hold(M::QuiverModuli)
  bounds = teleman_bounds(M.Q, M.d, M.theta)
  weights = all_weights_endomorphisms_universal_bundle(M.Q, M.d, M.theta)
  return all(maximum(weights[hn]) < bounds[hn] for hn in collect(keys(bounds)))
end

#####################################################################################
# Below are methods to work with Bundle objects that pertain to Teleman quantization
#####################################################################################

"""
    set_teleman_weights!(F::Oscar.AbstractBundle, weights::Dict{HNType,Vector{Int}})

Set the Teleman weights of the bundle `F` to the given dictionary.

This is used to assign Teleman weights to IntersectionTheory bundles via Oscar's
attribute system.
"""
function set_teleman_weights!(F::Oscar.AbstractBundle, weights::Dict{HNType,Vector{Int}})
  r = Int(Oscar.rank(F))
  !all(length(v) == r for v in values(weights)) &&
    throw(ArgumentError("Weights are not consistent with rank."))
  Oscar.set_attribute!(F, :teleman_weights, weights)
  return F
end
