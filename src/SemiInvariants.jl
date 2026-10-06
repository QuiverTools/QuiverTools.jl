# Cauchy decomposition assigns a partition to each arrow. At a vertex, the
# Schur products of outgoing and incoming partitions must differ by det^weight.

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

function _si_vertex_coefficient(
  context::_SchurContext, outgoing, incoming, weight::Int
)
  rank = context.rank
  rank == 0 && return big(1)
  sum(_si_size, outgoing; init=0) - sum(_si_size, incoming; init=0) ==
  weight * rank || return big(0)
  rank == 1 && return big(1)
  isempty(incoming) && return _si_rectangular_coefficient(context, outgoing, weight)
  isempty(outgoing) && return _si_rectangular_coefficient(context, incoming, -weight)

  # A single factor on one side fixes the partition needed on the other.
  if length(incoming) == 1 || length(outgoing) == 1
    fixed, factors, shift = if length(incoming) == 1
      (incoming[1], outgoing, weight)
    else
      (outgoing[1], incoming, -weight)
    end
    target = [_si_part(fixed, i) + shift for i in 1:rank]
    any(<(0), target) && return big(0)
    return _si_schur_coefficient(context, factors, _si_trim(target))
  end

  outgoing_product = _si_schur_product(context, outgoing)
  incoming_product = _si_schur_product(context, incoming)
  result = big(0)
  for (lambda, coefficient) in outgoing_product
    shifted = [_si_part(lambda, i) - weight for i in 1:rank]
    any(<(0), shifted) && continue
    result += coefficient * get(incoming_product, _si_trim(shifted), big(0))
  end
  return result
end

mutable struct _SIProblem
  dimensions::Vector{Int}
  weight::Vector{Int}
  order::Vector{Int}
  edges::Vector{Tuple{Int,Int}}
  outgoing::Vector{Vector{Int}}
  incoming::Vector{Vector{Int}}
  assigned::Vector{SIPartition}
  schur::Dict{Int,_SchurContext}
  vertex_cache::Dict{Tuple{Tuple,Tuple,Int,Int},BigInt}
  total::BigInt
end

function _SIProblem(Q::Quiver, dimensions::Vector{Int}, weight::Vector{Int})
  order = _si_topological_order(Q)
  edges = [(i, j) for (i, j) in arrows(Q) if dimensions[i] > 0 && dimensions[j] > 0]
  outgoing = [Int[] for _ in dimensions]
  incoming = [Int[] for _ in dimensions]
  for (index, (source, target)) in enumerate(edges)
    push!(outgoing[source], index)
    push!(incoming[target], index)
  end
  return _SIProblem(
    dimensions, weight, order, edges, outgoing, incoming,
    Vector{SIPartition}(undef, length(edges)), Dict{Int,_SchurContext}(),
    Dict{Tuple{Tuple,Tuple,Int,Int},BigInt}(), big(0),
  )
end

function _si_local_coefficient(problem::_SIProblem, vertex::Int)
  outgoing = Tuple(problem.assigned[edge] for edge in problem.outgoing[vertex])
  incoming = Tuple(problem.assigned[edge] for edge in problem.incoming[vertex])
  rank, weight = problem.dimensions[vertex], problem.weight[vertex]
  return get!(problem.vertex_cache, (outgoing, incoming, rank, weight)) do
    context = get!(problem.schur, rank) do
      _SchurContext(rank)
    end
    _si_vertex_coefficient(context, outgoing, incoming, weight)
  end
end

function _si_assign_outgoing!(
  problem::_SIProblem, position::Int, edge_position::Int,
  remaining::Int, coefficient::BigInt,
)
  vertex = problem.order[position]
  outgoing = problem.outgoing[vertex]
  if edge_position > length(outgoing)
    if remaining == 0
      local_coefficient = _si_local_coefficient(problem, vertex)
      local_coefficient == 0 ||
        _si_visit!(problem, position + 1, coefficient * local_coefficient)
    end
    return nothing
  end
  edge = outgoing[edge_position]
  target = problem.edges[edge][2]
  rank = min(problem.dimensions[vertex], problem.dimensions[target])
  context = get!(problem.schur, rank) do
    _SchurContext(rank)
  end
  for size in 0:remaining
    for partition in _si_partitions(context, size)
      problem.assigned[edge] = partition
      _si_assign_outgoing!(problem, position, edge_position + 1,
        remaining - size, coefficient)
    end
  end
  return nothing
end

# A determinant at a source fixes one outgoing partition, or makes two of
# them complementary in a rectangle. These are general Schur identities.
function _si_source!(problem::_SIProblem, position::Int, coefficient::BigInt)
  vertex = problem.order[position]
  isempty(problem.incoming[vertex]) && problem.dimensions[vertex] > 0 &&
  problem.weight[vertex] >= 0 || return false
  outgoing = problem.outgoing[vertex]
  rank, width = problem.dimensions[vertex], problem.weight[vertex]
  if length(outgoing) == 1
    edge = only(outgoing)
    target = problem.edges[edge][2]
    if width == 0 || rank <= problem.dimensions[target]
      problem.assigned[edge] = width == 0 ? () : Tuple(fill(width, rank))
      _si_visit!(problem, position + 1, coefficient)
    end
    return true
  end
  if length(outgoing) == 2
    first, second = outgoing
    first_rank = min(rank, problem.dimensions[problem.edges[first][2]])
    second_rank = min(rank, problem.dimensions[problem.edges[second][2]])
    _si_each_rectangle_partition(rank, width) do part
      lambda = _si_trim(copy(part))
      mu = _si_complement(lambda, width, rank)
      length(lambda) <= first_rank && length(mu) <= second_rank || return nothing
      problem.assigned[first], problem.assigned[second] = lambda, mu
      _si_visit!(problem, position + 1, coefficient)
    end
    return true
  end
  return false
end

function _si_visit!(problem::_SIProblem, position::Int, coefficient::BigInt)
  if position > length(problem.order)
    problem.total += coefficient
    return nothing
  end
  vertex = problem.order[position]
  required =
    problem.weight[vertex] * problem.dimensions[vertex] +
    sum(edge -> _si_size(problem.assigned[edge]), problem.incoming[vertex]; init=0)
  required < 0 && return nothing
  _si_source!(problem, position, coefficient) && return nothing
  _si_assign_outgoing!(problem, position, 1, required, coefficient)
  return nothing
end

function _si_generic_dimension(Q::Quiver, dimensions::Vector{Int}, weight::Vector{Int})
  problem = _SIProblem(Q, dimensions, weight)
  _si_visit!(problem, 1, big(1))
  return problem.total
end

function _si_reciprocal_problem(Q::Quiver, dimensions::Vector{Int}, weight::Vector{Int})
  alpha = zeros(Int, length(dimensions))
  for vertex in _si_topological_order(Q)
    alpha[vertex] =
      weight[vertex] +
      sum(Q.adjacency[source, vertex] * alpha[source]
          for source in eachindex(dimensions))
  end
  all(>=(0), alpha) || return nothing
  dual_weight = -euler_matrix(Q) * dimensions
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
  all(>=(0), d) || throw(ArgumentError("dimension vectors must be nonnegative"))
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
