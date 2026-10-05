# Dimensions of polynomial semi-invariant weight spaces for acyclic quivers.
# For an arrow i -> j, coordinate functions transform as V_i ⊗ V_j^*.

const SIPartition = Tuple{Vararg{Int}}

_part(p::SIPartition, i::Int) = i <= length(p) ? p[i] : 0
_partition_size(p::SIPartition) = sum(p; init=0)

function _si_trim(parts::Vector{Int})
  while !isempty(parts) && last(parts) == 0
    pop!(parts)
  end
  return Tuple(parts)
end

function _si_partitions(n::Int, rank::Int, cache)
  return get!(cache, (n, rank)) do
    n == 0 && return SIPartition[()]
    return SIPartition[Tuple(p) for p in partitions(n) if length(p) <= rank]
  end
end

function _si_topological_order(Q::Quiver)
  n = n_vertices(Q)
  indegrees = [indegree(Q, i) for i in 1:n]
  ready = [i for i in 1:n if indegrees[i] == 0]
  order = Int[]
  while !isempty(ready)
    i = pop!(ready)
    push!(order, i)
    for j in 1:n
      indegrees[j] -= Q.adjacency[i, j]
      if Q.adjacency[i, j] > 0 && indegrees[j] == 0
        push!(ready, j)
      end
    end
  end
  length(order) == n ||
    throw(ArgumentError("semi-invariant weight spaces require an acyclic quiver"))
  return order
end

function _si_horizontal_strip(outer::SIPartition, inner::SIPartition)
  for i in 1:max(length(outer), length(inner))
    _part(outer, i) >= _part(inner, i) >= _part(outer, i + 1) || return false
  end
  return true
end

function _si_each_rectangle_partition(f, height::Int, width::Int)
  part = zeros(Int, height)
  function visit(i::Int, upper::Int)
    if i > height
      f(part)
      return nothing
    end
    for value in 0:upper
      part[i] = value
      visit(i + 1, value)
    end
  end
  visit(1, width)
end

function _si_rectangle_square_shape(width::Int, part::Vector{Int})
  height = length(part)
  return _si_trim(
    vcat(
      [width + part[i] for i in 1:height],
      [width - part[i] for i in height:-1:1],
    ),
  )
end

function _si_rectangle_square(parts, rank::Int)
  length(parts) == 2 && parts[1] == parts[2] || return nothing
  rect = parts[1]
  !isempty(rect) && all(x == rect[1] for x in rect) || return nothing
  result = Dict{SIPartition,BigInt}()
  _si_each_rectangle_partition(length(rect), rect[1]) do part
    shape = _si_rectangle_square_shape(rect[1], part)
    length(shape) <= rank && (result[shape] = big(1))
  end
  return result
end

# The LR tableau is read from right to left in each row, starting at the top.
function _si_lr_coefficient(
  lambda::SIPartition, mu::SIPartition, nu::SIPartition, rank::Int, cache
)
  key = (lambda, mu, nu, rank)
  return get!(cache, key) do
    _partition_size(lambda) + _partition_size(mu) == _partition_size(nu) || return big(0)
    length(nu) <= rank || return big(0)
    all(_part(nu, i) >= _part(lambda, i) for i in 1:rank) || return big(0)
    isempty(mu) && return lambda == nu ? big(1) : big(0)
    isempty(lambda) && return mu == nu ? big(1) : big(0)
    if length(mu) == 1
      return _si_horizontal_strip(nu, lambda) ? big(1) : big(0)
    end
    if length(lambda) == 1
      return _si_horizontal_strip(nu, mu) ? big(1) : big(0)
    end
    if length(lambda) == rank && all(x == lambda[1] for x in lambda)
      return all(_part(nu, i) == lambda[1] + _part(mu, i) for i in 1:rank) ?
             big(1) : big(0)
    end
    if length(mu) == rank && all(x == mu[1] for x in mu)
      return all(_part(nu, i) == mu[1] + _part(lambda, i) for i in 1:rank) ?
             big(1) : big(0)
    end

    cells = Tuple{Int,Int}[]
    for i in 1:length(nu)
      for j in nu[i]:-1:(_part(lambda, i) + 1)
        push!(cells, (i, j))
      end
    end
    values = Dict{Tuple{Int,Int},Int}()
    used = zeros(Int, length(mu))
    function fill_cell(pos::Int)::BigInt
      pos > length(cells) && return big(1)
      row, col = cells[pos]
      right = get(values, (row, col + 1), length(mu))
      above = get(values, (row - 1, col), 0)
      result = big(0)
      for entry in (above + 1):right
        entry <= length(mu) || continue
        used[entry] < mu[entry] || continue
        used[entry] += 1
        if entry == 1 || used[entry - 1] >= used[entry]
          values[(row, col)] = entry
          result += fill_cell(pos + 1)
          delete!(values, (row, col))
        end
        used[entry] -= 1
      end
      return result
    end
    return fill_cell(1)
  end
end

function _si_schur_product(parts, rank::Int, partition_cache, lr_cache)
  square = _si_rectangle_square(parts, rank)
  square === nothing || return square
  pair, paired_product = nothing, nothing
  for i in 1:length(parts)
    for j in (i + 1):length(parts)
      candidate = _si_rectangle_square((parts[i], parts[j]), rank)
      if candidate !== nothing
        pair, paired_product = (i, j), candidate
        break
      end
    end
    pair === nothing || break
  end
  result = pair === nothing ? Dict{SIPartition,BigInt}(() => big(1)) :
           paired_product
  for (index, part) in enumerate(parts)
    pair !== nothing && index in pair && continue
    next = Dict{SIPartition,BigInt}()
    for (lambda, coefficient) in result
      degree = _partition_size(lambda) + _partition_size(part)
      for nu in _si_partitions(degree, rank, partition_cache)
        all(_part(nu, i) >= max(_part(lambda, i), _part(part, i)) for i in 1:rank) ||
          continue
        c = _si_lr_coefficient(lambda, part, nu, rank, lr_cache)
        c == 0 && continue
        next[nu] = get(next, nu, big(0)) + coefficient * c
      end
    end
    result = next
  end
  return result
end

function _si_schur_coefficient(parts, target::SIPartition, rank::Int,
  partition_cache, lr_cache)
  length(parts) == 0 && return isempty(target) ? big(1) : big(0)
  length(parts) == 1 && return parts[1] == target ? big(1) : big(0)
  length(parts) == 2 && return _si_lr_coefficient(parts[1], parts[2], target,
    rank, lr_cache)
  return get(_si_schur_product(parts, rank, partition_cache, lr_cache), target, big(0))
end

function _si_complement(part::SIPartition, width::Int, rank::Int)
  length(part) <= rank && _part(part, 1) <= width || return nothing
  return _si_trim([width - _part(part, rank + 1 - i) for i in 1:rank])
end

function _si_rectangular_coefficient(
  parts, rank::Int, width::Int, partition_cache, lr_cache
)
  width < 0 && return big(0)
  target = width == 0 ? () : Tuple(fill(width, rank))
  length(parts) == 0 && return width == 0 ? big(1) : big(0)
  length(parts) == 1 && return parts[1] == target ? big(1) : big(0)
  if length(parts) == 2
    return _si_complement(parts[1], width, rank) == parts[2] ? big(1) : big(0)
  end
  if length(parts) == 3
    complement = _si_complement(parts[3], width, rank)
    complement === nothing && return big(0)
    return _si_lr_coefficient(parts[1], parts[2], complement, rank, lr_cache)
  end
  row_index = findfirst(part -> length(part) == 1, parts)
  factors = if row_index === nothing
    parts
  else
    [parts[index] for index in eachindex(parts) if index != row_index]
  end
  middle = length(factors) ÷ 2
  left = _si_schur_product(factors[1:middle], rank, partition_cache, lr_cache)
  right = _si_schur_product(factors[(middle + 1):end], rank, partition_cache, lr_cache)
  result = big(0)
  for (lambda, a) in left
    complement = _si_complement(lambda, width, rank)
    complement === nothing && continue
    if row_index === nothing
      result += a * get(right, complement, big(0))
    else
      row_size = _partition_size(parts[row_index])
      for (mu, b) in right
        _partition_size(complement) - _partition_size(mu) == row_size || continue
        _si_horizontal_strip(complement, mu) && (result += a * b)
      end
    end
  end
  return result
end

function _si_vertex_coefficient(
  outgoing, incoming, rank::Int, weight::Int, partition_cache, lr_cache, vertex_cache
)
  rank == 0 && return big(1)
  sum(_partition_size, outgoing; init=0) - sum(_partition_size, incoming; init=0) ==
  weight * rank || return big(0)
  rank == 1 && return big(1)
  key = (Tuple(outgoing), Tuple(incoming), rank, weight)
  return get!(vertex_cache, key) do
    isempty(incoming) && return _si_rectangular_coefficient(
      outgoing, rank, weight, partition_cache, lr_cache
    )
    isempty(outgoing) && return _si_rectangular_coefficient(
      incoming, rank, -weight, partition_cache, lr_cache
    )
    if length(incoming) == 1
      target = [_part(incoming[1], i) + weight for i in 1:rank]
      any(x < 0 for x in target) && return big(0)
      return _si_schur_coefficient(outgoing, _si_trim(target), rank,
        partition_cache, lr_cache)
    end
    if length(outgoing) == 1
      target = [_part(outgoing[1], i) - weight for i in 1:rank]
      any(x < 0 for x in target) && return big(0)
      return _si_schur_coefficient(incoming, _si_trim(target), rank,
        partition_cache, lr_cache)
    end
    out = _si_schur_product(outgoing, rank, partition_cache, lr_cache)
    inn = _si_schur_product(incoming, rank, partition_cache, lr_cache)
    result = big(0)
    for (lambda, a) in out
      shifted = [_part(lambda, i) - weight for i in 1:rank]
      any(x < 0 for x in shifted) && continue
      result += a * get(inn, _si_trim(shifted), big(0))
    end
    return result
  end
end

function _si_generic_dimension(Q::Quiver, d::Vector{Int}, weight::Vector{Int})
  order = _si_topological_order(Q)
  edges = [(i, j) for (i, j) in arrows(Q) if d[i] > 0 && d[j] > 0]
  out_edges = [Int[] for _ in 1:length(d)]
  in_edges = [Int[] for _ in 1:length(d)]
  for (index, (i, j)) in enumerate(edges)
    push!(out_edges[i], index)
    push!(in_edges[j], index)
  end
  assigned = Vector{SIPartition}(undef, length(edges))
  partition_cache = Dict{Tuple{Int,Int},Vector{SIPartition}}()
  lr_cache = Dict{Any,BigInt}()
  vertex_cache = Dict{Any,BigInt}()
  total = big(0)

  function visit_vertex(position::Int, coefficient::BigInt)
    if position > length(order)
      total += coefficient
      return nothing
    end
    v = order[position]
    outgoing = out_edges[v]
    incoming = in_edges[v]
    required =
      weight[v] * d[v] + sum(
        index -> _partition_size(assigned[index]), incoming; init=0
      )
    required < 0 && return nothing
    # A determinant at a pure source with one arrow fixes its partition.
    # With two arrows, the partitions are complements in a rectangle.
    if isempty(incoming) && d[v] > 0 && weight[v] >= 0
      if length(outgoing) == 1
        edge = only(outgoing)
        target = edges[edge][2]
        weight[v] == 0 || d[v] <= d[target] || return nothing
        assigned[edge] = weight[v] == 0 ? () : Tuple(fill(weight[v], d[v]))
        visit_vertex(position + 1, coefficient)
        return nothing
      elseif length(outgoing) == 2
        first, second = outgoing
        first_rank = min(d[v], d[edges[first][2]])
        second_rank = min(d[v], d[edges[second][2]])
        _si_each_rectangle_partition(d[v], weight[v]) do part
          lambda = _si_trim(copy(part))
          mu = _si_complement(lambda, weight[v], d[v])
          length(lambda) <= first_rank && length(mu) <= second_rank || return nothing
          assigned[first], assigned[second] = lambda, mu
          visit_vertex(position + 1, coefficient)
        end
        return nothing
      end
    end
    function assign_outgoing(index::Int, remaining::Int)
      if index > length(outgoing)
        remaining == 0 || return nothing
        factor = _si_vertex_coefficient(
          [assigned[e] for e in outgoing], [assigned[e] for e in incoming],
          d[v], weight[v], partition_cache, lr_cache, vertex_cache,
        )
        factor == 0 || visit_vertex(position + 1, coefficient * factor)
        return nothing
      end
      local edge_index = outgoing[index]
      local target_vertex = edges[edge_index][2]
      for size in 0:remaining
        for part in _si_partitions(size, min(d[v], d[target_vertex]), partition_cache)
          assigned[edge_index] = part
          assign_outgoing(index + 1, remaining - size)
        end
      end
    end
    assign_outgoing(1, required)
  end
  visit_vertex(1, big(1))
  return total
end

function _si_reciprocal_problem(Q::Quiver, d::Vector{Int}, weight::Vector{Int})
  alpha = zeros(Int, length(d))
  for v in _si_topological_order(Q)
    alpha[v] = weight[v] + sum(Q.adjacency[u, v] * alpha[u] for u in 1:length(d))
  end
  all(x >= 0 for x in alpha) || return nothing
  dual_weight = -euler_matrix(Q) * d
  return alpha, Vector{Int}(dual_weight)
end

"""
    semi_invariant_dimension(Q::Quiver, d, weight; method=:lr)

Return the dimension of the polynomial semi-invariant weight space
``\\mathrm{SI}(Q,d)_{\\mathrm{weight}}`` as a `BigInt`.
The quiver must be acyclic. A coordinate on an arrow ``i\\to j`` has
determinant weight positive at ``i`` and negative at ``j``.
Weights at zero-dimensional vertices are ignored.

The default `method=:lr` uses the Cauchy and Littlewood--Richardson rules.
Rectangular Schur squares, rectangle complements, and Pieri's rule are applied
whenever their partition shapes occur, independently of the quiver.
Set `method=:reciprocity` to use Derksen--Weyman reciprocity
(MR1758751, Corollary 1). Reciprocity can be slower when its new dimension
vector is larger. This function concerns polynomial semi-invariants;
identifying them with sections of a line bundle on a quiver moduli space
requires the corresponding ample-stability hypothesis.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> semi_invariant_dimension(Q, [1, 1], [2, -2])
6

julia> semi_invariant_dimension(Q, [1, 1], [2, -2]; method=:reciprocity)
6

julia> semi_invariant_dimension(Quiver([0 1 1; 0 0 1; 0 0 0]), [1, 1, 1], [2, 0, -2])
3
```
"""
function semi_invariant_dimension(
  Q::Quiver, d::AbstractVector{<:Integer}, weight::AbstractVector{<:Integer};
  method::Symbol=:lr,
)
  n = n_vertices(Q)
  length(d) == n && length(weight) == n ||
    throw(DimensionMismatch("dimension and weight vectors must have one entry per vertex"))
  all(x >= 0 for x in d) || throw(ArgumentError("dimension vectors must be nonnegative"))
  method in (:lr, :reciprocity) ||
    throw(ArgumentError("method must be :lr or :reciprocity"))
  dimensions = Int.(d)
  sigma = [dimensions[i] == 0 ? 0 : Int(weight[i]) for i in 1:n]
  _si_topological_order(Q)
  sum(sigma[i] * dimensions[i] for i in 1:n) == 0 || return big(0)
  if method == :reciprocity
    reciprocal = _si_reciprocal_problem(Q, dimensions, sigma)
    reciprocal === nothing && return big(0)
    alpha, dual_weight = reciprocal
    return _si_generic_dimension(Q, alpha, dual_weight)
  end
  return _si_generic_dimension(Q, dimensions, sigma)
end
