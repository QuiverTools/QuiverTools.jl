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

# TODO is it arbitrary? is it not the Bourbaki orientation?
"""
    dynkin_quiver(type, n)

Construct the Dynkin quiver, with arbitrary orientation of the arrows.

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
    if n == 6
      return Quiver(
        [
          0 1 0 0 0 0
          0 0 1 0 0 0
          0 0 0 1 1 0
          0 0 0 0 0 0
          0 0 0 0 0 1
          0 0 0 0 0 0
        ],
        "Dynkin quiver of type E6",
      )
    elseif n == 7
      return Quiver(
        [
          0 1 0 0 0 0 0
          0 0 1 0 0 0 0
          0 0 0 1 1 0 0
          0 0 0 0 0 0 0
          0 0 0 0 0 1 0
          0 0 0 0 0 0 1
          0 0 0 0 0 0 0
        ],
        "Dynkin quiver of type E7",
      )
    elseif n == 8
      return Quiver(
        [
          0 1 0 0 0 0 0 0
          0 0 1 0 0 0 0 0
          0 0 0 1 1 0 0 0
          0 0 0 0 0 0 0 0
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

      )
    elseif n == 8
      return Quiver(
        [
          0 1 0 0 0 0 0 0 0
          0 0 1 0 0 0 0 0 0
          0 0 0 1 1 0 0 0 0
          0 0 0 0 0 0 0 0 0
          0 0 0 0 0 1 0 0 0
          0 0 0 0 0 0 1 0 0
          0 0 0 0 0 0 0 1 0
          0 0 0 0 0 0 0 0 1
          0 0 0 0 0 0 0 0 0
        ],
        "Dynkin quiver of type E8",
      )
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

""""
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
