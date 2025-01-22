##############
# Constructors
##############

export kronecker_quiver,
  loop_quiver,
  subspace_quiver,
  three_vertex_quiver,
  cyclic_quiver,
  bipartite_quiver,
  opposite_quiver,
  double_quiver,
  dynkin_quiver

export kronecker_moduli

"""
    kronecker_quiver(m::Int)

Construct the Kronecker quiver with `m` vertices.

# Input

- `m`: (Default = 2) The number of arrows in the Kronecker quiver.

# Output

The Kronecker quiver with `m` vertices.

# Examples

```jldoctest
julia> kronecker_quiver(3)
3-Kronecker quiver, with adjacency matrix [0 3; 0 0]
```
"""
function kronecker_quiver(m::Int=2)
  return Quiver([0 m; 0 0], string(m) * "-Kronecker quiver")
end

"""
    three_vertex_quiver(m12::Int, m13::Int, m23::Int)

Construct the three-vertex quiver with the given edge counts.

# Input

- `m12`: The number of arrows from vertex 1 to vertex 2.
- `m13`: The number of arrows from vertex 1 to vertex 3.
- `m23`: The number of arrows from vertex 2 to vertex 3.

# Output

A three-vertex quiver with the specified arrows.

# Examples

```jldoctest
julia> three_vertex_quiver(1, 2, 3)
Acyclic 3-vertex quiver, with adjacency matrix [0 1 2; 0 0 3; 0 0 0]
```
"""
function three_vertex_quiver(m12::Int, m13::Int, m23::Int)
  return Quiver([0 m12 m13; 0 0 m23; 0 0 0], "Acyclic 3-vertex quiver")
end

"""
    loop_quiver(m::Int)

Construct the loop quiver with `m` vertices.

# Input

- `m`: The number of vertices in the loop quiver.

# Output

The loop quiver with `m` vertices.

# Examples

```jldoctest
julia> loop_quiver(4)
4-loop quiver, with adjacency matrix [4;;]
```
"""
function loop_quiver(m::Int)
  return Quiver(Matrix{Int}(reshape([m], 1, 1)), string(m) * "-loop quiver")
end

"""
    subspace_quiver(m::Int)

Construct the subspace quiver with `m` vertices.

# Input

- `m`: The number of subspace-vertices.

# Output

The subspace quiver with `m` subspaces.

# Examples

```jldoctest
julia> subspace_quiver(3)
3-subspace quiver, with adjacency matrix [0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 0 0]
```
"""
function subspace_quiver(m::Int)
  A = zeros(Int, m + 1, m + 1)
  for i in 1:m
    A[i, m + 1] = 1
  end
  return Quiver(A, string(m) * "-subspace quiver")
end

function dynkin_quiver(Tn::String)
  type = Tn[1:(end - 1)]
  n = parse(Int, Tn[end])
  return dynkin_quiver(type, n)
end

# TODO is it arbitrary? is it not the Bourbaki orientation?
"""
    dynkin_quiver(type, n)

Construct the Dynkin quiver, with arbitrary orientation of the arrows.

# Examples

```jldoctest
julia> dynkin_quiver("D", 4)
Dynkin quiver of type D4, with adjacency matrix [0 1 0 0; 0 0 1 1; 0 0 0 0; 0 0 0 0]
```
"""
function dynkin_quiver(type::String, n::Int)
  if type == "A"
    if !(n >= 1)
      throw(ArgumentError("$n is out of bounds for type $type."))
    end
    if n == 1
      return Quiver([[1]], "Dynkin quiver of type A1")
    else
      M = zeros(Int, n, n)
      for i in 1:(n - 1)
        M[i, i + 1] = 1
      end
      return Quiver(M, "Dynkin quiver of type A$n")
    end
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
          0 1 0 0 0 0 0
          0 0 1 0 0 0 0
          0 0 0 1 1 0 0
          0 0 0 0 0 0 0
          0 0 0 0 0 1 0
          0 0 0 0 0 0 1
          0 0 0 0 0 0 0
        ],
        "Dynkin quiver of type E6",
      )
    elseif n == 7
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
        "Dynkin quiver of type E7",
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
cyclic quiver on 4 vertices, with adjacency matrix [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0]
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
julia> bipartite_quiver(2, 3).adjacency
5×5 StaticArraysCore.SMatrix{5, 5, Int64, 25} with indices SOneTo(5)×SOneTo(5):
 0  0  1  1  1
 0  0  1  1  1
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
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
2-Kronecker quiver, with adjacency matrix [0 2; 0 0]

julia> opposite_quiver(Q)
opposite of 2-Kronecker quiver, with adjacency matrix [0 0; 2 0]
```
"""
opposite_quiver(Q::Quiver) =
  Quiver(Matrix{Int}(transpose(Q.adjacency)), "opposite of " * Q.name)

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
double of 2-Kronecker quiver, with adjacency matrix [0 2; 2 0]
```
"""
double_quiver(Q::Quiver) =
  Quiver(Q.adjacency + Matrix{Int}(transpose(Q.adjacency)), "double of " * Q.name)

"""
    kronecker_moduli(m::Int, d::Int, e::Int)

Construct the Kronecker moduli space with `m` vertices and dimension vector `(d, e)`.

# Examples

```jldoctest
julia> kronecker_moduli(3, 2, 3)
Quiver moduli space defined as follows:
 - quiver: 3-Kronecker quiver, with adjacency matrix [0 3; 0 0]
 - dimension vector: [2, 3]
 - stability parameter: [9, -6]
 - condition: semistable
```
"""
function kronecker_moduli(m::Int, d::Int, e::Int)
  return QuiverModuliSpace(kronecker_quiver(m), [d, e])
end
