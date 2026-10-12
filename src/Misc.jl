######
# Misc
######

# Shared validator for public quiver-setting functions that accept `(Q, d)` directly.
function __check_dimension_vector(Q::Quiver, d::AbstractVector{Int})
  length(d) == n_vertices(Q) ||
    throw(ArgumentError("dimension vector must have length $(n_vertices(Q))"))
  all(>=(0), d) || throw(ArgumentError("dimension vector must be non-negative"))
  return nothing
end

"""
    identity_matrix(n::Int)

Return the identity matrix of size `n`.
"""
@memoize function identity_matrix(n::Int)
  out = zeros(Int, n, n)
  for i in 1:n
    out[i, i] = 1
  end
  return coerce_matrix(out)
end

"""
    diagonal(m::AbstractMatrix{Int})

Return the diagonal matrix with the diagonal of `m` as its diagonal.
"""
function diagonal(m::AbstractMatrix{Int})
  n = size(m, 1)
  return map(ind -> ind[1] == ind[2] ? m[ind...] : 0, Iterators.product(1:n, 1:n))
end

"""
    diagonal(v::AbstractVector)

Return a square matrix with diagonal `v`.
"""
function diagonal(v::AbstractVector)
  n = length(v)
  return map(ind -> ind[1] == ind[2] ? v[ind[1]] : 0, Iterators.product(1:n, 1:n))
end

function support(d::AbstractVector{Int})
  return findall(x -> x != 0, d)
end

"""
    support_subquiver(Q::Quiver, d::AbstractVector{Int})

Restrict `(Q, d)` to the support of `d`: the full subquiver on the vertices where `d` is
nonzero, together with `d` restricted to those vertices.

The moduli space `M(Q, d)` is isomorphic to `M(support_subquiver(Q, d)...)`, since a zero
entry `dᵢ = 0` forces `Vᵢ = 0` in every representation. Results that only hold for a
full-support dimension vector (e.g. the Mukai index, see
[[MR4352662](https://mathscinet.ams.org/mathscinet-getitem?mr=4352662)]) must therefore be
computed on this restriction rather than on `(Q, d)` directly.

For internal use only.
"""
function support_subquiver(Q::Quiver, d::AbstractVector{Int})
  supp = support(d)
  return Quiver(Matrix(Q.adjacency[supp, supp])), d[supp]
end

########################################################################################
# Technical tools
########################################################################################

"""
    zero_vector(n::Int)

Create a zero vector of length `n`.

# Examples

```jldoctest
julia> QuiverTools.zero_vector(3) == [0, 0, 0]
true
```
"""
zero_vector(n::Int) = coerce_vector(zeros(Int, n))

"""
  zero_vector(Q::Quiver)

Create the zero dimension vector for the quiver `Q`.

# Examples

```jldoctest
julia> QuiverTools.zero_vector(kronecker_quiver(3)) == [0, 0]
true
```
"""
zero_vector(Q::Quiver) = zero_vector(n_vertices(Q))

"""
    thin_dimension_vector(Q::Quiver)

Create the thin dimension vector for a given quiver `Q`.

# Examples:

```jldoctest
julia> Q = kronecker_quiver(3);

julia> QuiverTools.thin_dimension_vector(Q) == [1, 1]
true
```
"""
thin_dimension_vector(Q::Quiver) = coerce_vector(ones(Int, n_vertices(Q)))

"""
    all_subdimension_vectors(d::AbstractVector{Int}; nonzero::Bool=false, strict::Bool=false)

Compute all subdimension vectors of a given dimension vector `d`.

# Input

- `d::AbstractVector{Int}`: The input dimension vector.
- `nonzero::Bool=false`: whether to exclude the zero vector.
- `strict::Bool=false`: whether to exclude the input vector `d`.

# Output

- An array of all subdimension vectors of `d`, with or without the zero vector and `d`.

# Examples

```jldoctest
julia> QuiverTools.all_subdimension_vectors([2, 3])
12-element Vector{Vector{Int64}}:
 [0, 0]
 [1, 0]
 [2, 0]
 [0, 1]
 [1, 1]
 [2, 1]
 [0, 2]
 [1, 2]
 [2, 2]
 [0, 3]
 [1, 3]
 [2, 3]

julia> QuiverTools.all_subdimension_vectors([2, 3]; nonzero=true)
11-element Vector{Vector{Int64}}:
 [1, 0]
 [2, 0]
 [0, 1]
 [1, 1]
 [2, 1]
 [0, 2]
 [1, 2]
 [2, 2]
 [0, 3]
 [1, 3]
 [2, 3]

julia> QuiverTools.all_subdimension_vectors([2, 3]; nonzero=true, strict=true)
10-element Vector{Vector{Int64}}:
 [1, 0]
 [2, 0]
 [0, 1]
 [1, 1]
 [2, 1]
 [0, 2]
 [1, 2]
 [2, 2]
 [0, 3]
 [1, 3]

julia> QuiverTools.all_subdimension_vectors([0, 0, 0]; nonzero=true, strict=true)
Vector{Int64}[]
```
"""
function all_subdimension_vectors(
  d::AbstractVector{Int};
  nonzero::Bool=false,
  strict::Bool=false,
)
  @assert length(d) > 0 "Input vector must have positive length."
  n = length(d)
  if all(di == 0 for di in d)
    return (nonzero || strict) ? Vector{Int}[] : [zeros(Int, n)]
  end
  total = prod(di + 1 for di in d)
  out = Vector{Vector{Int}}(undef, total)
  digits = zeros(Int, n)
  for idx in 1:total
    out[idx] = copy(digits)
    # increment as a mixed-radix counter, coordinate 1 varying fastest, which
    # reproduces the order of the previous recursive `vcat` builder
    for k in 1:n
      digits[k] += 1
      if digits[k] <= d[k]
        break
      else
        digits[k] = 0
      end
    end
  end
  nonzero && popfirst!(out)  # drop the zero vector (first element)
  strict && pop!(out)        # drop d itself (last element)
  return out
end

"""
    is_subdimension_vector(e::AbstractVector{Int}, d::AbstractVector{Int})

Check whether vector `e` is a subdimension of vector `d`.

# Examples

```jldoctest
julia> QuiverTools.is_subdimension_vector([1, 1], [2, 3])
true

julia> QuiverTools.is_subdimension_vector([1, 1], [1, 1])
true

julia> QuiverTools.is_subdimension_vector([1, 2], [1, 1])
false
```
"""
function is_subdimension_vector(e::AbstractVector{Int}, d::AbstractVector{Int})
  return all(ei <= di for (ei, di) in zip(e, d))
end

"""
    unit_vector(n::Int, i::Int)

Return a vector of length `n` with a `1` at index `i` and `0` elsewhere.

# Arguments

- `n::Int`: The length of the unit vector.
- `i::Int`: The index at which to place the `1` in the unit vector.

# Examples

```jldoctest
julia> QuiverTools.unit_vector(3, 2) == [0, 1, 0]
true
```
"""
function unit_vector(n::Int, i::Int)
  v = zeros(Int, n)
  v[i] = 1
  return coerce_vector(v)
end

"""
    unit_vector(Q::Quiver, i::Int)

Return a dimension vector for the quiver `Q` with a `1` at index `i` and `0` elsewhere.

# Arguments

- `Q::Quiver`: The input quiver.
- `i::Int`: The index at which to place the `1` in the unit vector.

# Examples
```jldoctest
julia> Q = kronecker_quiver(3);

julia> QuiverTools.unit_vector(Q, 2) == [0, 1]
true
```
"""
unit_vector(Q::Quiver, i::Int) = unit_vector(n_vertices(Q), i)

coerce_vector(v) = v
coerce_matrix(m) = SMatrix{size(m, 1),size(m, 1)}(m)
