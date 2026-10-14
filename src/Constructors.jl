##############
# Constructors
##############

"""
    kronecker_quiver(m::Int)

Construct the Kronecker quiver with `m` arrows.

# Input

- `m`: (Default = 2) The number of arrows.

# Output

The Kronecker quiver with `m` arrows.

# Examples

```jldoctest
julia> kronecker_quiver(3)
3-Kronecker quiver
```
"""
function kronecker_quiver(m::Int=2)
  return Quiver([0 m; 0 0], string(m) * "-Kronecker quiver")
end

"""
    three_vertex_quiver(m12::Int, m13::Int, m23::Int)

Construct the three-vertex quiver with the given arrow counts.

# Input

- `m12`: The number of arrows from vertex 1 to vertex 2.
- `m13`: The number of arrows from vertex 1 to vertex 3.
- `m23`: The number of arrows from vertex 2 to vertex 3.

# Output

A three-vertex quiver with the specified arrows.

# Examples

```jldoctest
julia> three_vertex_quiver(1, 2, 3)
Acyclic 3-vertex quiver
```
"""
function three_vertex_quiver(m12::Int, m13::Int, m23::Int)
  return Quiver([0 m12 m13; 0 0 m23; 0 0 0], "Acyclic 3-vertex quiver")
end

"""
    loop_quiver(m::Int)

Construct the loop quiver with `m` arrows.

# Input

- `m`: The number of arrows.

# Output

The loop quiver with `m` arrows.

# Examples

```jldoctest
julia> loop_quiver(4)
4-loop quiver
```
"""
function loop_quiver(m::Int)
  return Quiver(reshape([m], 1, 1), string(m) * "-loop quiver")
end

"""
    jordan_quiver(m::Int=1)

Construct the Jordan quiver, i.e. the quiver with one vertex and `m` loops.

This is a synonym of `loop_quiver`: for `m = 1` it is the Jordan quiver,
otherwise the generalized Jordan quiver with `m` loops.

# Input

- `m`: (Default = 1) The number of loops.

# Output

The quiver with one vertex and `m` loops.

# Examples

```jldoctest
julia> jordan_quiver()
Jordan quiver

julia> jordan_quiver(3)
generalized Jordan quiver with 3 loops
```
"""
function jordan_quiver(m::Int=1)
  name = m == 1 ? "Jordan quiver" : "generalized Jordan quiver with $m loops"
  return Quiver(reshape([m], 1, 1), name)
end

"""
    subspace_quiver(m::Int)

Construct the subspace quiver with `m + 1` vertices.

# Input

- `m`: The number of subspace-vertices.

# Output

The subspace quiver with `m` subspaces.

# Examples

```jldoctest
julia> subspace_quiver(3)
3-subspace quiver

julia> n_vertices(subspace_quiver(3))
4
```
"""
function subspace_quiver(m::Int)
  A = zeros(Int, m + 1, m + 1)
  for i in 1:m
    A[i, m + 1] = 1
  end
  return Quiver(A, string(m) * "-subspace quiver")
end

"""
    star_quiver(m::Int)

Synonym for [`subspace_quiver`](@ref).
"""
star_quiver(m::Int) = subspace_quiver(m)

"""
    generalized_subspace_quiver(m::Int, K::AbstractVector{Int})

Construct the generalized subspace quiver with `m + 1` vertices and `K[i]` arrows
from the `i`-th source to the sink.

# Input

- `m`: The number of subspace-vertices (sources).
- `K`: A vector of length `m`; `K[i]` is the number of arrows from source `i` to the sink.

# Output

The generalized subspace quiver with `m` sources and multiplicities `K`.

# Examples

```jldoctest
julia> generalized_subspace_quiver(3, [1, 2, 3])
a generalized 3-subspace quiver
```
"""
function generalized_subspace_quiver(m::Int, K::AbstractVector{Int})
  length(K) == m || throw(ArgumentError("K must have length m = $m"))
  A = zeros(Int, m + 1, m + 1)
  for i in 1:m
    A[i, m + 1] = K[i]
  end
  return Quiver(A, "a generalized $m-subspace quiver")
end

"""
    thickened_subspace_quiver(m::Int, k::Int)

Construct the thickened subspace quiver with `m + 1` vertices and `k` arrows
from each of the `m` sources to the sink.

# Input

- `m`: The number of subspace-vertices (sources).
- `k`: The number of arrows from each source to the sink.

# Output

The thickened subspace quiver with `m` sources, each with `k` arrows to the sink.

# Examples

```jldoctest
julia> thickened_subspace_quiver(3, 2)
thickened subspace quiver with 3 sources and multiplicity 2
```
"""
function thickened_subspace_quiver(m::Int, k::Int)
  A = generalized_subspace_quiver(m, fill(k, m)).adjacency
  return Quiver(A, "thickened subspace quiver with $m sources and multiplicity $k")
end

"""
    framed_quiver(Q::Quiver, n::AbstractVector{Int})

Construct the framed quiver ``\\widehat{Q}`` of `Q` with framing datum `n`.

A new framing vertex ``i_0`` is prepended as the **first** vertex, together with `n[i]`
arrows ``i_0 \\to i`` for every vertex `i` of `Q`. See the framing construction in
[arXiv:2607.12895](https://arxiv.org/abs/2607.12895).

# Input

- `Q`: a quiver.
- `n`: a vector of length `n_vertices(Q)`; `n[i]` is the number of framing arrows to vertex `i`.

# Output

The framed quiver, with the framing vertex as vertex `1`.

# Examples

```jldoctest
julia> framed_quiver(kronecker_quiver(2), [0, 1])
framing of 2-Kronecker quiver
```
"""
function framed_quiver(Q::Quiver, n::AbstractVector{Int})
  length(n) == n_vertices(Q) ||
    throw(ArgumentError("length of n must equal the number of vertices"))
  N = n_vertices(Q)
  A = zeros(Int, N + 1, N + 1)
  A[2:end, 2:end] .= Q.adjacency
  A[1, 2:end] .= n
  return Quiver(A, "framing of " * Q.name)
end

"""
    coframed_quiver(Q::Quiver, n::AbstractVector{Int})

Construct the coframed quiver of `Q` with coframing datum `n`.

A new coframing vertex ``i_0`` is appended as the **last** vertex, together with `n[i]`
arrows ``i \\to i_0`` for every vertex `i` of `Q`. This is the linear dual of
[`framed_quiver`](@ref); see [arXiv:2607.12895](https://arxiv.org/abs/2607.12895).

# Input

- `Q`: a quiver.
- `n`: a vector of length `n_vertices(Q)`; `n[i]` is the number of coframing arrows from vertex `i`.

# Output

The coframed quiver, with the coframing vertex as the last vertex.

# Examples

```jldoctest
julia> coframed_quiver(kronecker_quiver(2), [0, 1])
coframing of 2-Kronecker quiver
```
"""
function coframed_quiver(Q::Quiver, n::AbstractVector{Int})
  length(n) == n_vertices(Q) ||
    throw(ArgumentError("length of n must equal the number of vertices"))
  N = n_vertices(Q)
  A = zeros(Int, N + 1, N + 1)
  A[1:N, 1:N] .= Q.adjacency
  A[1:N, N + 1] .= n
  return Quiver(A, "coframing of " * Q.name)
end

# Split a Dynkin label like "A3" or "D10" into its letter type and integer rank.
function _parse_dynkin_label(Tn::String)
  m = match(r"^([A-Za-z]+)([0-9]+)$", Tn)
  if isnothing(m)
    throw(ArgumentError("$Tn is not a valid Dynkin label, e.g. \"A3\" or \"D10\"."))
  end
  return String(m.captures[1]), parse(Int, m.captures[2])
end

"""
    dynkin_quiver(Tn::String)

Construct the Dynkin quiver from a type string such as `"A3"` or `"D10"`.

See `dynkin_quiver(type, n)` for details.

# Examples

```jldoctest
julia> dynkin_quiver("D10")
Dynkin quiver of type D10
```
"""
function dynkin_quiver(Tn::String)
  type, n = _parse_dynkin_label(Tn)
  return dynkin_quiver(type, n)
end

"""
    dynkin_quiver(type, n)

Construct the Dynkin quiver on the Bourbaki-numbered vertices,
oriented lexicographically (arrows go from lower- to higher-numbered vertices).

Supported types are `"A"` (`n ≥ 1`), `"D"` (`n ≥ 3`) and `"E"` (`n ∈ {6, 7, 8}`).

# Examples

```jldoctest
julia> dynkin_quiver("D", 4)
Dynkin quiver of type D4

julia> dynkin_quiver("A", 1)
Dynkin quiver of type A1

julia> n_vertices(dynkin_quiver("E", 6))
6
```
"""
function dynkin_quiver(type::String, n::Int)
  if type == "A"
    if !(n >= 1)
      throw(ArgumentError("$n is out of bounds for type $type."))
    end
    M = zeros(Int, n, n)
    for i in 1:(n - 1)
      M[i, i + 1] = 1
    end
    return Quiver(M, "Dynkin quiver of type A$n")
  elseif type == "D"
    if !(n >= 3)
      throw(ArgumentError("$n is out of bounds for type $type."))
    end
    M = zeros(Int, n, n)
    for i in 1:(n - 2)
      M[i, i + 1] = 1
    end
    M[n - 2, n] = 1

    return Quiver(M, "Dynkin quiver of type D$n")
  elseif type == "E"
    if !(n in [6, 7, 8])
      throw(ArgumentError("$n is out of bounds for type $type."))
    end
    # Bourbaki numbering: the chain is 1—3—4—⋯—n, with vertex 2 attached to vertex 4
    if n == 6
      return Quiver(
        [
          0 0 1 0 0 0
          0 0 0 1 0 0
          0 0 0 1 0 0
          0 0 0 0 1 0
          0 0 0 0 0 1
          0 0 0 0 0 0
        ],
        "Dynkin quiver of type E6",
      )
    elseif n == 7
      return Quiver(
        [
          0 0 1 0 0 0 0
          0 0 0 1 0 0 0
          0 0 0 1 0 0 0
          0 0 0 0 1 0 0
          0 0 0 0 0 1 0
          0 0 0 0 0 0 1
          0 0 0 0 0 0 0
        ],
        "Dynkin quiver of type E7",
      )
    elseif n == 8
      return Quiver(
        [
          0 0 1 0 0 0 0 0
          0 0 0 1 0 0 0 0
          0 0 0 1 0 0 0 0
          0 0 0 0 1 0 0 0
          0 0 0 0 0 1 0 0
          0 0 0 0 0 0 1 0
          0 0 0 0 0 0 0 1
          0 0 0 0 0 0 0 0
        ],
        "Dynkin quiver of type E8",
      )
    end
  else
    throw(ArgumentError("$type is not a valid ADE Dynkin type."))
  end
end

"""
    extended_dynkin_quiver(Tn::String)

Construct the extended (affine) Dynkin quiver from a type string such as `"A3"` or `"D10"`.

See `extended_dynkin_quiver(type, n)` for details.

# Examples

```jldoctest
julia> extended_dynkin_quiver("D10")
Extended Dynkin quiver of type D10
```
"""
function extended_dynkin_quiver(Tn::String)
  type, n = _parse_dynkin_label(Tn)
  return extended_dynkin_quiver(type, n)
end

"""
    extended_dynkin_quiver(type, n)

Construct the extended (affine) Dynkin quiver of type `type` with `n + 1` vertices,
oriented lexicographically (arrows go from lower- to higher-numbered vertices).

The vertices follow the Bourbaki numbering of the affine diagram, shifted by one
so that the affine node ``0`` becomes vertex `1` and the Bourbaki node ``k``
becomes vertex `k + 1`.

Supported types are `"A"` (`n ≥ 1`), `"D"` (`n ≥ 4`) and `"E"` (`n ∈ {6, 7, 8}`).

# Examples

```jldoctest
julia> extended_dynkin_quiver("A", 1) == kronecker_quiver()
true

julia> extended_dynkin_quiver("A", 2) == three_vertex_quiver(1, 1, 1)
true

julia> extended_dynkin_quiver("D", 4)
Extended Dynkin quiver of type D4

julia> n_vertices(extended_dynkin_quiver("E", 6))
7
```
"""
function extended_dynkin_quiver(type::String, n::Int)
  if type == "A"
    n >= 1 || throw(ArgumentError("$n is out of bounds for type $type."))
    if n == 1
      return Quiver([0 2; 0 0], "Extended Dynkin quiver of type A1")
    end
    M = zeros(Int, n + 1, n + 1)
    for i in 1:n
      M[i, i + 1] = 1
    end
    M[1, n + 1] = 1
    return Quiver(M, "Extended Dynkin quiver of type A$n")
  elseif type == "D"
    n >= 4 || throw(ArgumentError("$n is out of bounds for type $type."))
    M = zeros(Int, n + 1, n + 1)
    M[1, 3] = 1
    M[2, 3] = 1
    for i in 3:(n - 2)
      M[i, i + 1] = 1
    end
    M[n - 1, n] = 1
    M[n - 1, n + 1] = 1
    return Quiver(M, "Extended Dynkin quiver of type D$n")
  elseif type == "E"
    # Bourbaki numbering of the finite diagram shifted by one (vertex k + 1 is
    # the Bourbaki node k), with the affine node as vertex 1
    if n == 6
      return Quiver(
        [
          0 0 1 0 0 0 0
          0 0 0 1 0 0 0
          0 0 0 0 1 0 0
          0 0 0 0 1 0 0
          0 0 0 0 0 1 0
          0 0 0 0 0 0 1
          0 0 0 0 0 0 0
        ],
        "Extended Dynkin quiver of type E6",
      )
    elseif n == 7
      return Quiver(
        [
          0 1 0 0 0 0 0 0
          0 0 0 1 0 0 0 0
          0 0 0 0 1 0 0 0
          0 0 0 0 1 0 0 0
          0 0 0 0 0 1 0 0
          0 0 0 0 0 0 1 0
          0 0 0 0 0 0 0 1
          0 0 0 0 0 0 0 0
        ],
        "Extended Dynkin quiver of type E7",
      )
    elseif n == 8
      return Quiver(
        [
          0 0 0 0 0 0 0 0 1
          0 0 0 1 0 0 0 0 0
          0 0 0 0 1 0 0 0 0
          0 0 0 0 1 0 0 0 0
          0 0 0 0 0 1 0 0 0
          0 0 0 0 0 0 1 0 0
          0 0 0 0 0 0 0 1 0
          0 0 0 0 0 0 0 0 1
          0 0 0 0 0 0 0 0 0
        ],
        "Extended Dynkin quiver of type E8",
      )
    else
      throw(ArgumentError("$n is out of bounds for type $type."))
    end
  else
    throw(ArgumentError("$type is not a valid ADE Dynkin type."))
  end
end

"""
    cyclic_quiver(n)

Construct the cyclic quiver on `n` vertices.

# Examples

```jldoctest
julia> cyclic_quiver(4)
cyclic quiver on 4 vertices
```
"""
function cyclic_quiver(n::Int)
  if n < 1
    throw(ArgumentError("$n must be greater than 0"))
  end
  A = zeros(Int, n, n)
  for i in 1:(n - 1)
    A[i, i + 1] = 1
  end
  A[n, 1] = 1
  return Quiver(A, "cyclic quiver on $n vertices")
end

"""
    bipartite_quiver(m, n)

Construct the bipartite quiver on `m` and `n` vertices.

# Input

- `m`: The number of vertices in the first part.
- `n`: The number of vertices in the second part.

# Output

The bipartite quiver with `m + n` vertices.

# Examples

```jldoctest
julia> print(bipartite_quiver(2, 3).adjacency)
[0 0 1 1 1; 0 0 1 1 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
```
"""
function bipartite_quiver(m::Int, n::Int)
  if m < 1 || n < 1
    throw(ArgumentError("m and n must be greater than 0"))
  end
  A = zeros(Int, m + n, m + n)
  for i in 1:m
    for j in (m + 1):(m + n)
      A[i, j] = 1
    end
  end
  return Quiver(A, "bipartite quiver on $m and $n vertices")
end

"""
    opposite_quiver(Q::Quiver)

Construct the opposite quiver.

The opposite quiver has the same vertices,
and an arrow ``j \\to i`` for every arrow  ``i \\to j`` in the original quiver.

# Input

- `Q::Quiver`: The quiver to be reversed.

# Output

The quiver with the same vertices and reversed arrows.

# Examples

```jldoctest
julia> Q = kronecker_quiver()
2-Kronecker quiver

julia> opposite_quiver(Q)
opposite of 2-Kronecker quiver
```
"""
opposite_quiver(Q::Quiver) = Quiver(
  transpose(Q.adjacency), "opposite of " * Q.name
)

"""
    double_quiver(Q::Quiver)

Construct the double of a quiver.

The double of a quiver has the same vertices,
and for every arrow ``i \\to j`` in the original quiver, an arrow ``j\\to i`` is added.
In other words, the double of a quiver is the quiver obtained by adding the transpose
of the adjacency matrix.

# Examples

```jldoctest
julia> Q = kronecker_quiver();

julia> double_quiver(Q)
double of 2-Kronecker quiver
```
"""
double_quiver(Q::Quiver) = Quiver(
  Q.adjacency + transpose(Q.adjacency), "double of " * Q.name
)

"""
    disjoint_union(Q1::Quiver, Q2::Quiver)

Construct the disjoint union of two quivers.

The disjoint union has the vertices and arrows of both quivers and no arrows between
them; its adjacency matrix is the block diagonal of the two adjacency matrices.

# Examples

```jldoctest
julia> disjoint_union(kronecker_quiver(3), loop_quiver(2))
disjoint union of 3-Kronecker quiver and 2-loop quiver
```
"""
function disjoint_union(Q1::Quiver, Q2::Quiver)
  n1, n2 = n_vertices(Q1), n_vertices(Q2)
  A = zeros(Int, n1 + n2, n1 + n2)
  A[1:n1, 1:n1] = Q1.adjacency
  A[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)] = Q2.adjacency
  name = if (!isempty(Q1.name) && !isempty(Q2.name))
    "disjoint union of $(Q1.name) and $(Q2.name)"
  else
    ""
  end
  return Quiver(A, name)
end

"""
    kronecker_moduli(m::Int, d::Int, e::Int)

Construct the Kronecker moduli space with `m` vertices and dimension vector `(d, e)`.

# Examples

```jldoctest
julia> kronecker_moduli(3, 2, 3)
Quiver moduli space defined as follows:
 - quiver: 3-Kronecker quiver
 - dimension vector: [2, 3]
 - stability parameter: [9, -6]
 - condition: semistable
```
"""
kronecker_moduli(m::Int, d::Int, e::Int) = QuiverModuliSpace(kronecker_quiver(m), [d, e])

"""
    subspace_quiver_moduli(m::Int, d::Int)

Construct the subspace quiver moduli space for `m` points on ``\\mathbb{P}^{d-1}``

# Examples

```jldoctest
julia> subspace_quiver_moduli(5, 2)
Quiver moduli space defined as follows:
 - quiver: 5-subspace quiver
 - dimension vector: [1, 1, 1, 1, 1, 2]
 - stability parameter: [2, 2, 2, 2, 2, -5]
 - condition: semistable
```
"""
subspace_quiver_moduli(m::Int, d::Int) = QuiverModuliSpace(
  subspace_quiver(m), vcat(repeat([1], m), [d])
)
