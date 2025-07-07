######
# Misc
######

"""
    identity_matrix(n::Int)

Return the identity matrix of size `n`.
"""
@memoize Dict identity_matrix(n::Int) = map(
  ind -> ind[1] == ind[2] ? 1 : 0, Iterators.product(1:n, 1:n)
)

"""
    diagonal(m::AbstractMatrix{Int})

Return the diagonal matrix with the diagonal of `m` as its diagonal.
"""
function diagonal(m::AbstractMatrix{Int})
  n = size(m)[1]
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
@memoize Dict zero_vector(n::Int) = coerce_vector(zeros(Int, n))

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

julia> QuiverTools.all_subdimension_vectors([0, 0, 0]; nonzero=true, strict=true)
StaticArraysCore.SVector{3, Int64}[]
```
"""
@memoize Dict function all_subdimension_vectors(
  d::AbstractVector{Int};
  nonzero::Bool=false,
  strict::Bool=false,
) #TODO should this be memoized at all?
  subdims = reshape(
    coerce_vector.(collect(Iterators.product(map(di -> 0:di, d)...))),
    prod(di + 1 for di in d),
  )
  all(di == 0 for di in d) && (nonzero || strict) && return deleteat!(subdims, 1)
  nonzero && deleteat!(subdims, 1)
  strict && pop!(subdims)
  return subdims
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
@memoize Dict function unit_vector(n::Int, i::Int)
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

coerce_vector(v::AbstractVector) = SVector{length(v)}(v)
coerce_vector(v::Tuple) = SVector{length(v)}(v)
coerce_vector(v::SVector) = v

coerce_matrix(m::AbstractMatrix) = SMatrix{size(m)...}(m)
coerce_matrix(m::SMatrix) = m
