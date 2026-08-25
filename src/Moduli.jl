####################################
# Methods to deal with quiver moduli
####################################

"""
    is_nonempty(M::QuiverModuli)

Checks if the quiver moduli is nonempty.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

# Output

- whether the moduli space is nonempty.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> is_nonempty(M)
true

julia> M = QuiverModuliSpace(Q, [2, 3], [-3, 2]);

julia> is_nonempty(M)
false
```
"""
function is_nonempty(M::QuiverModuli)
  if M.condition == "stable"
    return has_stables(M.Q, M.d, M.theta)
  elseif M.condition == "semistable"
    return has_semistables(M.Q, M.d, M.theta)
  end
end

"""
    is_coprime(M::QuiverModuli)

Checks if the stability parameter is coprime with the dimension vector,
i.e., if for all subdimension vectors ``e`` of ``d``, ``\\theta\\cdot e \\neq 0``.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

# Output

- whether the dimension vector `M.d` is theta-coprime for `M.theta`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> is_coprime(M)
true
```
"""
function is_coprime(M::QuiverModuli)
  return is_coprime(M.d, M.theta)
end

"""
    all_hn_types(M::QuiverModuli; unstable::Bool=false, ordered::Bool=true)

Returns all Harder-Narasimhan types of the moduli space.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

Keyword arguments:

- `unstable::Bool`: if `true`, returns only Harder-Narasimhan types
corresponding to unstable representations. Default is `false`
- `ordered::Bool`: if `true`, returns the Harder-Narasimhan types in
the order introduced by [[MR1974891](https://doi.org/10.1007/s00222-002-0273-4)].
Default is `true`.

# Output

- a list of Harder-Narasimhan types for the dimension vector and slope of `M`.

# Examples

The HN types for a 3-Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> all_hn_types(M)
8-element Vector{HNType}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> all_hn_types(M; unstable = true)
7-element Vector{HNType}:
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]
```
"""
function all_hn_types(M::QuiverModuli; unstable::Bool=false, ordered::Bool=true)
  return all_hn_types(M.Q, M.d, M.theta, M.denom; unstable=unstable, ordered=ordered)
end

"""
    is_hn_type(M::QuiverModuli, hn_type::HNType)

Checks if the given sequence of dimension vectors is a valid HN type for
the moduli space.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `hn_type::HNType`: a Harder--Narasimhan type.

# Output

- whether the given sequence is a valid Harder-Narasimhan type for `M`.

# Examples

Some HN types for the 3-Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_hn_type(M, [[2, 3]])
true

julia> is_hn_type(M, [[1, 1], [1, 2]])
true

julia> is_hn_type(M, [[1, 2], [1, 1]])
false
```
"""
is_hn_type(M::QuiverModuli, hn_type::HNType) = is_hn_type(
  M.Q, M.d, hn_type, M.theta, M.denom
)
is_hn_type(M::QuiverModuli, hn_type::Vector{<:AbstractVector{Int}}) = is_hn_type(
  M, HNType(hn_type)
)

"""
    codimension_hn_stratum(M::QuiverModuli, hn_type::HNType)

Computes the codimension of the Harder-Narasimhan stratum
corresponding to the given HN type.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `hn_type::HNType`: a Harder--Narasimhan type

# Output

- the codimension of the Harder-Narasimhan stratum corresponding to the given HN type.

# Examples

Codimensions for the 3-Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> codimension_hn_stratum(M, [[2, 3]])
0

julia> codimension_hn_stratum(M, [[1, 1], [1, 2]])
3
```
"""
codimension_hn_stratum(M::QuiverModuli, hn_type::HNType) = codimension_hn_stratum(
  M.Q, hn_type
)
codimension_hn_stratum(M::QuiverModuli, hn_type::Vector{<:AbstractVector{Int}}) =
  codimension_hn_stratum(
    M, HNType(hn_type)
  )

"""
    codimension_unstable_locus(M::QuiverModuli)

Computes the codimension of the unstable locus in the parameter space.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

# Output

- the codimension of the unstable locus in the parameter space.

# Examples

Codimensions for the 3-Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> codimension_unstable_locus(M)
3
```

If the unstable locus is empty, the codimension is `Inf`:
```jldoctest
julia> Q = kronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 0]);

julia> all_hn_types(M; unstable=true)
HNType[]

julia> codimension_unstable_locus(M)
Inf
```
"""
function codimension_unstable_locus(M::QuiverModuli)
  hn_types = all_hn_types(M; unstable=true)
  isempty(hn_types) && return Inf
  return minimum(codimension_hn_stratum(M, hn_type) for hn_type in hn_types)
end

"""
    all_luna_types(M::QuiverModuli; stable::Bool = true)

Returns all Luna types of the moduli space.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

Keyword arguments:

- `stable::Bool`: if `false`, excludes the stable Luna type. Default is `true`

# Output

- a list of Luna types for the dimension vector and slope of `M`. Each entry is a
  [`LunaType`](@ref): a dictionary whose keys are the distinct dimension vectors
  ``\\mathbf{d}^k`` occurring in the type and whose values are lists of the
  multiplicities with which they occur. For example `Dict([1, 1] => [2, 1])` means the
  dimension vector `[1, 1]` appears twice, once with multiplicity `2` and once with
  multiplicity `1`. See [`LunaType`](@ref) for the full description.

# Examples

Luna types for a 3-Kronecker quiver:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [3, 3]);

julia> all_luna_types(M)
5-element Vector{LunaType}:
 Dict([3, 3] => [1])
 Dict([1, 1] => [1], [2, 2] => [1])
 Dict([1, 1] => [3])
 Dict([1, 1] => [1, 2])
 Dict([1, 1] => [1, 1, 1])
```
"""
function all_luna_types(M::QuiverModuli; stable::Bool=true)
  return all_luna_types(M.Q, M.d, M.theta, M.denom; stable=stable)
end

"""
    all_luna_types(Q::Quiver, d, theta, denom; stable=true)

Computes all the possible Luna types for the given data.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter. Defaults to `canonical_stability(Q, d)`.
- `denom::Function`: a function defining the denominator of the slope. Defaults to `sum`.

Keyword arguments:

- `stable::Bool`: if `false`, excludes the stable Luna type. Default is `true`

# Output

- a list of Luna types. Each entry is a [`LunaType`](@ref): a dictionary whose keys are
  the distinct dimension vectors ``\\mathbf{d}^k`` occurring in the type and whose values
  are lists of the multiplicities with which they occur. For example
  `Dict([1, 1] => [2, 1])` means the dimension vector `[1, 1]` appears twice, once with
  multiplicity `2` and once with multiplicity `1`. See [`LunaType`](@ref) for the full
  description.


# Examples

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [3, 3]);

julia> all_luna_types(M)
5-element Vector{LunaType}:
 Dict([3, 3] => [1])
 Dict([1, 1] => [1], [2, 2] => [1])
 Dict([1, 1] => [3])
 Dict([1, 1] => [1, 2])
 Dict([1, 1] => [1, 1, 1])

julia> X = QuiverModuliSpace(Q, [2, 3]);

julia> all_luna_types(X)
1-element Vector{LunaType}:
 Dict([2, 3] => [1])
```
"""
@memoize Dict function all_luna_types(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum;
  stable::Bool=true,
)

  # treat the zero case separately
  all(di == 0 for di in d) && return [LunaType(Dict(d => [1]))]

  # subdimensions with the same slope as d
  μ = slope(d, theta, denom)
  same_slope = all_subdimension_vectors(d; nonzero=true, strict=true)
  filter!(e -> slope(e, theta, denom) == μ, same_slope)
  filter!(e -> has_stables(Q, e, theta, denom), same_slope)

  luna_types = LunaType[]
  for e in same_slope
    for luna_type in all_luna_types(Q, d - e, theta, denom; stable=true)
      if haskey(luna_type, e)
        for i in eachindex(luna_type[e])
          push!(luna_types, __add_and_return(luna_type, e, i))
        end
        # A second distinct stable summand of dimension `e` only exists if `e` admits
        # more than one isomorphism class of stable representation, i.e. its stable
        # locus is positive-dimensional (`euler_form(Q, e, e) <= 0`). A rigid `e`
        # (`euler_form(Q, e, e) == 1`) has a unique stable representation, so repeating
        # it would produce a Luna type that no representation realizes (issue #24).
        euler_form(Q, e, e) <= 0 &&
          push!(luna_types, __add_and_return_noniso(luna_type, e))
      else
        push!(luna_types, __add_and_return_new(luna_type, e))
      end
    end
  end

  stable && has_stables(Q, d, theta, denom) &&
    pushfirst!(luna_types, LunaType(Dict(d => [1])))

  # sort the multiplicities vectors, so that we can apply `unique!` to remove duplicates
  #
  # Example: this avoids having both `Dict([1, 1] => [2, 1])` and `Dict([1, 1] => [1, 2])`
  map(
    lt -> begin
      for key in keys(lt)
        sort!(lt[key])
      end
    end,
    luna_types,
  )
  return unique!(luna_types)
end

"""
    __add_and_return(luna_type, e, i)

Returns a new Luna type obtained by increasing the multiplicity of the `i`-th copy
the subdimension vector `e` by 1, in the given Luna type.

Internal use only.
"""
function __add_and_return(luna_type, e, i)
  new_luna_type = deepcopy(luna_type)
  new_luna_type[e][i] += 1
  return new_luna_type
end

"""
    __add_and_return_noniso(luna_type, e)

Returns a new Luna type obtained by adding a new copy
of the subdimension vector `e` in the given Luna type.
This assumes that there exists at least one more representation
that is not isomorphic to the present one.

Internal use only.
"""
function __add_and_return_noniso(luna_type, e)
  new_luna_type = deepcopy(luna_type)
  pushfirst!(new_luna_type[e], 1)
  return new_luna_type
end

"""
    __add_and_return_new(luna_type, e)

Returns a new Luna type obtained by adding a subdimension vector `e`
with multiplicity 1 to the given Luna type.

Internal use only.
"""
function __add_and_return_new(luna_type, e)
  new_luna_type = deepcopy(luna_type)
  new_luna_type[e] = [1]
  return new_luna_type
end

"""
    is_luna_type(M::QuiverModuli, tau)

Checks if the given tau is a valid Luna type for `M`.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `tau::Dict{AbstractVector{Int}, Vector{Int}}`: a candidate Luna type for `M`, encoded
  as a dictionary of dimension vectors to lists of multiplicities; see [`LunaType`](@ref).

# Output

- whether the given tau is a valid Luna type for `M`.

# Examples

Nontrivial Luna types for the 3-Kronecker quiver:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [3, 3]);

julia> l = Dict([1, 1] => [1], [2, 2] => [1]);

julia> is_luna_type(M, l)
true
```

The zero dimensional case:
```jldoctest
julia> Q = kronecker_quiver(3); X = QuiverModuliSpace(Q, [0, 0]);

julia> is_luna_type(X, Dict([0, 0] => [1]))
true
```
"""
function is_luna_type(M::QuiverModuli, tau)
  if sum(M.d) == 0
    return tau == Dict(M.d => [1])
  end

  ks = collect(keys(tau))
  isempty(ks) && return false
  if !all(
    e -> length(e) == n_vertices(M.Q) && all(>=(0), e) && any(>(0), e),
    ks,
  )
    return false
  end
  if !all(e -> !isempty(tau[e]) && all(>(0), tau[e]), ks)
    return false
  end
  # each key `e` contributes `sum(tau[e])` copies of `e` (one per multiplicity in its list)
  if sum(sum(tau[e]) * e for e in ks) != M.d
    return false
  end
  if !all(slope(e, M.theta, M.denom) == slope(M.d, M.theta, M.denom) for e in ks)
    return false
  end

  if !all(has_stables(M.Q, e, M.theta, M.denom) for e in ks)
    return false
  end
  # A rigid stable representation is unique up to isomorphism, so its dimension
  # vector cannot encode several distinct stable summands in one Luna type.
  return all(e -> length(tau[e]) == 1 || euler_form(M.Q, e, e) <= 0, ks)
end

"""
    dimension_of_luna_stratum(M::QuiverModuli, tau)

Computes the dimension of the Luna stratum corresponding to the given Luna type in the
moduli space.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `tau::Dict{AbstractVector{Int}, Vector{Int}}`: a Luna type for `M`, encoded as a
  dictionary of dimension vectors to lists of multiplicities; see [`LunaType`](@ref).

# Output

- the dimension of the Luna stratum corresponding to the given Luna type.

# Examples

```jldoctest
julia> Q = kronecker_quiver(2); M = QuiverModuliSpace(Q, [2, 2], [1, -1]);

julia> luna = all_luna_types(M)
2-element Vector{LunaType}:
 Dict([1, 1] => [2])
 Dict([1, 1] => [1, 1])

julia> [dimension_of_luna_stratum(M, tau) for tau in luna]
2-element Vector{Int64}:
 1
 2

julia> M = QuiverModuliSpace(Q, [0, 0]);

julia> dimension_of_luna_stratum(M, Dict([0, 0] => [1]))
0
```
"""
function dimension_of_luna_stratum(M::QuiverModuli, tau)
  is_luna_type(M, tau) ||
    throw(DomainError(tau, "not a Luna type for the given moduli problem"))
  # the formula below would give 1 for the zero dimension vector
  sum(M.d) == 0 && return 0
  return sum(length(tau[e]) * (1 - euler_form(M.Q, e, e)) for e in collect(keys(tau)))
end

"""
    local_quiver_setting(M::QuiverModuli, tau)

Returns the local quiver and dimension vector for the given Luna type.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `tau::Dict{AbstractVector{Int}, Vector{Int}}`: a Luna type for `M`.

# Output

- a named tuple `(Q, d, summands)` containing the local quiver, its dimension vector,
  and the dimension vectors of the stable summands, one for each vertex of the local
  quiver, ordered compatibly with `d`.
"""
function local_quiver_setting(M::QuiverModuli, tau)
  if !is_luna_type(M, tau)
    throw(DomainError(tau, "not a Luna type for the given moduli problem"))
  end

  # one local vertex per distinct stable summand, i.e. per entry of each multiplicity list;
  # `summands` and `dloc` iterate the keys in the same order, so vertex `k` carries `dloc[k]`.
  summands = [e for e in keys(tau) for _m in tau[e]]
  s = length(summands)
  # the number of arrows from vertex k to vertex l is δ_{k,l} - ⟨d_k, d_l⟩, see MR1972892
  A = [
    (k == l ? 1 : 0) - euler_form(M.Q, summands[k], summands[l]) for k in 1:s, l in 1:s
  ]

  Qloc = Quiver(A)
  dloc = [m for e in keys(tau) for m in tau[e]]

  return (Q=Qloc, d=dloc, summands=summands)
end

# whether the local quiver setting of the Luna type is coregular, i.e., whether the
# moduli space is smooth along the corresponding stratum
function __is_smooth_stratum(M::QuiverModuli, tau)
  setting = local_quiver_setting(M, tau)
  return is_coregular(setting.Q, setting.d)
end

"""
    semistable_equals_stable(M::QuiverModuli)

Checks if stability and semistability are equivalent on the given moduli space.
In other words, checks if there are no properly semistable points in the representation
space.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

# Output

- whether every semistable representation is stable.

# Examples

If the dimension vector is coprime with the stability parameter, then semistability
and stability are equivalent:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> semistable_equals_stable(M)
true
```

However, this is not necessarily the case:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [3, 3]);

julia> semistable_equals_stable(M)
false
```
"""
function semistable_equals_stable(M::QuiverModuli)
  if is_coprime(M.d, M.theta) || !has_semistables(M.Q, M.d, M.theta, M.denom)
    return true
  end
  return length(all_luna_types(M; stable=false)) == 0
end

"""
    is_amply_stable(M::QuiverModuli)

Checks whether the dimension vector ``d`` is amply stable
with respect to the slope function `theta`/`denominator`.

This means that the codimension of the unstable locus
in the parameter space is at least ``2``.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

# Output

- `true` if the codimension of the unstable locus is at least `2`, `false` otherwise.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_amply_stable(M)
true
```
"""
function is_amply_stable(M::QuiverModuli)
  return codimension_unstable_locus(M) >= 2
end

"""
    is_strongly_amply_stable(M::QuiverModuli)

Checks whether the quiver moduli setup `M.Q, M.d` is strongly amply `M.theta`-stable,
see [[Definition 4.1, doi:10.5802/jep.312](https://doi.org/10.5802/jep.312)].


# Examples

Example 4.8 from [[doi:10.5802/jep.312](https://doi.org/10.5802/jep.312)]:
```jldoctest
julia> Q = Quiver("1-----2-3,1-3"); M = QuiverModuliSpace(Q, [4, 1, 4]);

julia> is_strongly_amply_stable(M)
false

julia> is_amply_stable(M)
true
```

Our favorite example of a strongly amply stable setup:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_strongly_amply_stable(M)
true
```
"""
function is_strongly_amply_stable(M::QuiverModuli)
  return is_strongly_amply_stable(M.Q, M.d, M.theta)
end

"""
    dimension(M::QuiverModuliStack)

Returns the dimension of the moduli stack.
This differs from the dimension of the moduli space by 1, as we do not quotient out
the stabilizer `` \\mathbb{G}``.

# Input

- `M::QuiverModuliStack`: a moduli stack of representations of a quiver.

# Output

- the dimension of the moduli stack.

# Examples

The dimension of the moduli stack of the 3-Kronecker quiver

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliStack(Q, [2, 3]);

julia> dimension(M)
5
```
"""
function dimension(M::QuiverModuliStack)
  if is_nonempty(M)
    return -euler_form(M.Q, M.d, M.d)
  end
  return -Inf
end

"""
    dimension(M::QuiverModuliSpace)

Returns the dimension of the moduli space.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the dimension of the moduli space.

# Examples

The dimension of the moduli space of the 3-Kronecker quiver:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> dimension(M)
6
```
"""
function dimension(M::QuiverModuliSpace)
  # dimension is independent of the linearization, so cache it on `M.chow`. In `unsafe`
  # mode this cache is pre-seeded with `1 - <d, d>`, so we never run the expensive
  # `has_stables` check (issue #20). Only a finite dimension is cached; the empty case
  # returns the `-Inf` sentinel, which is cheap to recompute and stays uncached.
  cached = M.chow._dimension
  cached !== nothing && return cached
  n = _dimension(M)
  n isa Int && setfield!(M.chow, :_dimension, n)
  return n
end

function _dimension(M::QuiverModuliSpace)
  # the zero representation is semistable, but not stable, for d = 0
  !is_connected(M.Q) &&
    throw(ArgumentError("Q is not connected, M has disjoint connected components."))

  if all(M.d .== 0)
    if M.condition == "semistable"
      return 0
    else
      return -Inf
    end
  end

  # if the stable locus is nonempty then the dimension is 1 - <d, d>
  if has_stables(M.Q, M.d, M.theta, M.denom)
    return 1 - euler_form(M.Q, M.d, M.d)
  end

  # if the stable locus is empty, the dimension is the maximum of the dimensions
  # of the Luna strata
  if M.condition == "stable"
    return -Inf
  elseif M.condition == "semistable"
    if has_semistables(M.Q, M.d, M.theta, M.denom)
      return maximum(
        dimension_of_luna_stratum(M, tau) for
        tau in all_luna_types(M.Q, M.d, M.theta, M.denom)
      )
    end
  end
  # the semistable locus is also empty
  return -Inf
end

"""
    is_smooth(M::QuiverModuliSpace)

Checks if the moduli space is smooth.

In the presence of properly semistable representations, the moduli space is
étale-locally isomorphic, around a polystable representation, to the affine quotient of
the corresponding local quiver setting near the zero representation, by
[[MR1972892](https://mathscinet.ams.org/mathscinet/relay-station?mr=1972892)].
Following the strategy of
[[Theorem 4.2, MR1929191](https://mathscinet.ams.org/mathscinet/relay-station?mr=1929191)],
the moduli space is thus smooth if and only if the local quiver setting of every Luna
type is coregular, which is checked using [`is_coregular`](@ref).

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- whether the moduli space is smooth.

# Examples

Setups with `d` `theta`-coprime are smooth:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_smooth(M)
true
```

For the 3-Kronecker quiver and `d = (3, 3)` the moduli space is singular, whereas for
`d = (2, 2)` and `d = (2, 4)` one gets ``\\mathbb{P}^5``, despite the presence of
properly semistable representations:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [3, 3]);

julia> is_smooth(M)
false

julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 2]);

julia> is_smooth(M)
true

julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 4]);

julia> is_smooth(M)
true
```
"""
function is_smooth(M::QuiverModuliSpace)
  if M.condition == "stable"
    return true
  elseif !has_properly_semistables(M.Q, M.d, M.theta, M.denom)
    return true
  end

  # smoothness at the polystable points of a Luna stratum is equivalent to
  # coregularity of its local quiver setting, by combining the étale-local description
  # of [MR1972892] with [Theorem 2.1, MR1929191]; this is the globalization of
  # [Theorem 4.2, MR1929191] to arbitrary stability parameters
  return all(tau -> __is_smooth_stratum(M, tau), all_luna_types(M))
end

"""
    codimension_singular_locus(M::QuiverModuliSpace)

Computes the codimension of the singular locus of the moduli space.

The singular locus is a union of Luna strata: all points of the stratum of a Luna type
are singular if the corresponding local quiver setting is not coregular, and smooth
otherwise, as in [`is_smooth`](@ref). Unlike for moduli of vector bundles on a curve,
the singular locus can be strictly smaller than the locus of properly semistable
representations, whose codimension is bounded by that of the singular locus.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the codimension of the singular locus, or `Inf` if the moduli space is smooth.

# Examples

The Segre cubic threefold, with its ten singular points:
```jldoctest
julia> M = QuiverModuliSpace(subspace_quiver(6), [1, 1, 1, 1, 1, 1, 2]);

julia> codimension_singular_locus(M)
3
```

For the 3-Kronecker quiver and `d = (2, 2)` the properly semistable locus is non-empty
yet the moduli space is smooth, whilst for `d = (3, 3)` there are singularities:
```jldoctest
julia> codimension_singular_locus(QuiverModuliSpace(kronecker_quiver(3), [2, 2]))
Inf

julia> codimension_singular_locus(QuiverModuliSpace(kronecker_quiver(3), [3, 3]))
3
```
"""
function codimension_singular_locus(M::QuiverModuliSpace)
  M.condition == "stable" && return Inf

  # the stratum of a Luna type consists of singular points if and only if its local
  # quiver setting is not coregular; the stable stratum is always smooth
  singular = filter(tau -> !__is_smooth_stratum(M, tau), all_luna_types(M))
  isempty(singular) && return Inf
  return dimension(M) - maximum(dimension_of_luna_stratum(M, tau) for tau in singular)
end

"""
    is_smooth(M::QuiverModuliStack)

Checks if the moduli stack is smooth.

This is always trus, as the quotient stack of a smooth variety is smooth.

# Input

- `M::QuiverModuliStack`: a moduli stack of representations of a quiver.

# Output

- true

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliStack(Q, [2, 3]);

julia> is_smooth(M)
true
```
"""
function is_smooth(M::QuiverModuliStack)
  return true
end

"""
    is_projective(M::QuiverModuli)

Checks if the moduli space is projective.

# Input

- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

# Output

- whether the moduli space is projective.

# Examples

The moduli space of the 3-Kronecker quiver is projective:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_projective(M)
true
```
"""
function is_projective(M::QuiverModuli)
  if is_acyclic(M.Q)
    M.condition == "semistable" && return true
    M.condition == "stable" && return !has_properly_semistables(M.Q, M.d, M.theta, M.denom)
  end

  SSP = semisimple_moduli_space(M)
  M.condition == "semistable" && return dimension(SSP) in [0, -Inf]
  M.condition == "stable" &&
    return (
      dimension(SSP) in [1, -Inf] && !has_properly_semistables(M.Q, M.d, M.theta, M.denom)
    )
end

"""
    semisimple_moduli_space(M::QuiverModuli)

Returns the moduli space with the zero stability parameter.

Any quiver moduli space is (quasi)projective-over-affine;
this is the affine base.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the moduli space with the zero stability parameter.

# Examples
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> dimension(semisimple_moduli_space(M))
0
```
"""
function semisimple_moduli_space(M::QuiverModuliSpace)
  return QuiverModuliSpace(M.Q, M.d, zero_vector(M.Q))
end
