######
# Misc
######

# this is the wheel reinvention department.
# I don't want to load the whole LinearAlgebra package just for this.
# TODO is there a good reason for not wanting to do this? it's a standard library package?
"""
    identity_matrix(n::Int)

Returns the identity matrix of size `n`.
"""
@memoize Dict identity_matrix(n::Int) =
  map(ind -> ind[1] == ind[2] ? 1 : 0, Iterators.product(1:n, 1:n))

"""
    diagonal(m::AbstractMatrix{Int})

Returns a copy of the input matrix `m` with the diagonal untouched,
and all other entries set to zero.
"""
function diagonal(m::AbstractMatrix{Int})
  n = size(m)[1]
  return map(ind -> ind[1] == ind[2] ? m[ind...] : 0, Iterators.product(1:n, 1:n))
end

"""
    diagonal(v::AbstractVector)

Returns a square matrix with diagonal `v`.
"""
function diagonal(v::AbstractVector)
  n = length(v)
  return map(ind -> ind[1] == ind[2] ? v[ind[1]] : 0, Iterators.product(1:n, 1:n))
end

########################################################################################
# Technical tools
########################################################################################

"""
    zero_vector(n::Int)

Create a zero vector of length `n`.

# Input

- `n::Int`: The length of the zero vector.

# Output

- A zero vector of length `n`.

EXAMPLE:

There is not much to it:
```jldoctest
julia> QuiverTools.zero_vector(3) == [0, 0, 0]
true
```
"""
@memoize Dict function zero_vector(n::Int)
  return coerce_vector(zeros(Int, n))
end

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

Compute the thin dimension vector for a given quiver `Q`.

# Input

- `Q::Quiver`: The input quiver.

# Output

- A vector of ones of length `n`.

# Examples:

There is not much to it:
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
12-element Vector{StaticArraysCore.SVector{2, Int64}}:
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
11-element Vector{StaticArraysCore.SVector{2, Int64}}:
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
10-element Vector{StaticArraysCore.SVector{2, Int64}}:
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
```
"""
@memoize Dict function all_subdimension_vectors(
  d::AbstractVector{Int};
  nonzero::Bool=false,
  strict::Bool=false,
)
  subdimension_vectors = coerce_vector.(collect(Iterators.product(map(di -> 0:di, d)...)))
  if nonzero
    subdimension_vectors = filter(e -> any(ei != 0 for ei in e), subdimension_vectors)
  end
  if strict
    subdimension_vectors = filter(e -> e != d, subdimension_vectors)
  end
  return filter(e -> true, subdimension_vectors)
end

"""
    is_subdimension_vector(e::AbstractVector{Int}, d::AbstractVector{Int})

Check if vector `e` is a subdimension of vector `d`.

# Input

- `e::AbstractVector{Int}` A vector of integers.
- `d::AbstractVector{Int}` A vector of integers.

# Output

whether `e` is a subdimension of `d`.

EXAMPLE:
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

# Input

- `n::Int`: The length of the unit vector.
- `i::Int`: The index at which to place the `1` in the unit vector.

# Output

A unit vector of length `n` with a `1` at index `i` and `0` elsewhere.

# Examples

```jldoctest
julia> QuiverTools.unit_vector(3, 2) == [0, 1, 0]
true
```
"""
@memoize Dict function unit_vector(n::Int, i::Int)
  v = zeros(Int, n)
  v[i] = 1
  return coerce_vector(v)
end

"""
    unit_vector(Q::Quiver, i::Int)

Return a dimension vector for the quiver `Q` with a `1` at index `i` and `0` elsewhere.

# Input

- `Q::Quiver`: The input quiver.
- `i::Int`: The index at which to place the `1` in the unit vector.

# Output

- A dimension vector for the quiver `Q` with a `1` at index `i` and `0` elsewhere.

# Examples
```jldoctest
julia> Q = kronecker_quiver(3);

julia> QuiverTools.unit_vector(Q, 2) == [0, 1]
true
```
"""
unit_vector(Q::Quiver, i::Int) = unit_vector(n_vertices(Q), i)

coerce_vector(v::AbstractVector) = SVector{length(v)}(v)
coerce_vector(v::Tuple) = SVector{length(v)}(v)
coerce_vector(v::SVector) = v

coerce_matrix(m::AbstractMatrix) = SMatrix{size(m)...}(m)
coerce_matrix(m::SMatrix) = m
