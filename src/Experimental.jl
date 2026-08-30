# Here goes code that is not meant to be in a release.
"""
    mutation(Q::Quiver, i::int; keep_name::Bool=false)


Perform the mutation of the quiver `Q` at vertex `i`. This is defined
in [[Definition 1.1, 2112.11812](https://arxiv.org/abs/2112.11812)].


# Input

- `Q::Quiver` a quiver.
- `i::Int` the index of the vertex at which to mutate.
- `keep_name::Bool` whether to keep the name of the quiver in the output. Default is `false`.

# Output

The mutated quiver.

# Examples

```jldoctest
julia> Q = Quiver("1---2,2---3,3---1");

julia> Q1 = mutation(Q, 3)
Quiver with adjacency matrix [0 0 3; 6 0 0; 0 3 0]

julia> Q2 = mutation(Q1, 2)
Quiver with adjacency matrix [0 6 0; 0 0 3; 15 0 0]

julia> K = kronecker_quiver(3); mutation(K, 1; keep_name=true)
mutation of 3-Kronecker quiver at vertex 1

julia> K = kronecker_quiver(3); K1 = mutation(K, 1; keep_name=true)
mutation of 3-Kronecker quiver at vertex 1

julia> K1.adjacency
2×2 Matrix{Int64}:
 0  0
 3  0
```
"""
function mutation(Q::Quiver, i::Int; keep_name::Bool=false)
  n = n_vertices(Q)
  @assert (i >= 1 && i <= n) "Index i must be between 1 and n_vertices(Q)"

  # entry (j, k) is # of paths j -> i -> k
  paths_through_i = map(jk -> Q.adjacency[jk[1], i] * Q.adjacency[i, jk[2]], Iterators.product(1:n, 1:n))

  # copy() returns an immutable object, we must built it anew
  new_adjacency = zeros(Int, size(Q.adjacency))
  map(ij -> new_adjacency[ij[1], ij[2]] = Q.adjacency[ij[1], ij[2]], Iterators.product(1:n, 1:n))

  for j in 1:n
    if j != i
      # reverse arrows between i and j
      new_adjacency[j, i] = Q.adjacency[i, j]
      new_adjacency[i, j] = Q.adjacency[j, i]
    end
  end

  # add paths through i
  for j in 1:n, k in 1:n
    if j != i && k != i
        new_adjacency[j, k] += paths_through_i[j, k]
    end
  end

  # remove 2-cycles
  for j in 1:(n-1)
    for k in (j+1):n
      if new_adjacency[j, k] > new_adjacency[k, j]
        new_adjacency[j, k] -= new_adjacency[k, j]
        new_adjacency[k, j] = 0
      else
        new_adjacency[k, j] -= new_adjacency[j, k]
        new_adjacency[j, k] = 0
      end
    end
  end

  !keep_name && return Quiver(new_adjacency)
  return Quiver(new_adjacency, "mutation of " * Q.name * " at vertex $i")
end


"""
    mutation(Q::Quiver, d::AbstractVector{Int}, i::Int; keep_name::Bool=false)

Perform the mutation of the quiver `Q` at vertex `i` and update the dimension vector `d` accordingly.
This is defined in [[Definition 1.1, 2112.11812](https://arxiv.org/abs/2112.11812)].

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector for `Q`.
- `i::Int` the index of the vertex at which to mutate.
- `keep_name::Bool` whether to keep the name of the quiver in the output. Default is `false`.

# Output

The mutated quiver and the updated dimension vector.

# Examples

```jldoctest
julia> Q = Quiver("1---2,2---3,3---1"); d = [1, 1, 1];

julia> Q1, d1 = mutation(Q, d, 3)
(Quiver with adjacency matrix [0 0 3; 6 0 0; 0 3 0], [1, 1, 2])
```
"""
function mutation(Q::Quiver, d::AbstractVector{Int}, i::Int; keep_name::Bool=false)
  new_Q = mutation(Q, i; keep_name=keep_name)
  new_d = copy(d)
  new_d[i] = maximum([sum(Q.adjacency[j, i] * d[j] for j in 1:n_vertices(Q)), sum(Q.adjacency[i, j] * d[j] for j in 1:n_vertices(Q))]) - d[i]
  return new_Q, new_d
end



function mori_mukai_invariant(M::QuiverModuliSpace)
  Kb_dual = dual(canonical_bundle(M))

  K = chern_class(Kb_dual)
  D = [QuiverTools.Oscar.gens(chow_ring(M))[i] for i in accumulate(+, [1; M.d])[1:end-1]];
  deleteat!(D, findfirst(x -> x != 0, M.chow.chi))

  return LinearAlgebraX.cofactor_det([div(-K * Di * Dj, point_class(M)) for Di in D, Dj in D])
end

function is_amply_stable_fast(M::QuiverModuliSpace)
  # we only need to enumerate HN types of length 2
  for e in QuiverTools.all_destabilizing_subdimension_vectors(M.d, M.theta)
    !has_semistables(M.Q, e, M.theta, M.denom) && continue
    !has_semistables(M.Q, M.d - e, M.theta, M.denom) && continue
    hn = HNType([e, M.d - e])
    QuiverTools.codimension_hn_stratum(M.Q, hn) == 1 && return false
  end
  return true
end



function framing(Q::Quiver; n::AbstractVector{Int}=ones(Int, n_vertices(Q)))
  @assert length(n) == n_vertices(Q) "Framing `n` must be of length `n_vertices(Q)`."
  A = zeros(Int, n_vertices(Q) + 1, n_vertices(Q) + 1)
  A[2:end, 2:end] = Q.adjacency
  for j in 1:n_vertices(Q)
    A[1, j+1] = n[j]
  end

  Q.name == "" && return Quiver(A)
  return Quiver(A, "framing of " * Q.name * " by " * string(n))
end

function framing(Q::Quiver, d::AbstractVector{Int}; n::AbstractVector{Int}= ones(Int, n_vertices(Q)))
  new_Q = framing(Q; n=n)
  new_d = vcat(1, d)
  return new_Q, new_d
end

"""
    general_stability(chamber)

Given a chamber of the VGIT fan,
return a stability parameter in its relative interior.

In practice this is a rescaled sum of the rays of the polyhedral cone input.
"""
function general_stability(chamber)
  out = sum(QuiverTools.Oscar.rays(chamber))

  c = lcm(denominator.(out)...)
  out .*= c

  l = gcd(out...)
  return Int.(out ./ l)
end




function framing_stability(Q, d, theta)
  throw(NotImplementedError("The stability parameter for smooth models needs to be implemented as in doi:10.1007/s00209-008-0401-y"))
end


"""
    is_QFV(Q::Quiver, d::AbstractVector{Int}, n::AbstractVector{Int})

Compute whether the framing of the datum `Q, d` by `n` is a quiver flag variety.
"""
function is_QFV(Q, d, n)
  F, df = framing(Q, d; n=n)

  nu = pushfirst!(- ones(Int, n_vertices(Q)), sum(d))
  theta_can = canonical_stability(F, df)
  return git_equivalent(F, df, nu, theta_can)
end

function find_framing(Q, d; bound::Int=10)
  # this function finds the smallest framing for which the smooth model is a quiver flag variety

  for e in QuiverTools.all_subdimension_vectors(bound*d; nonzero=true)
    !all(ei > 0 for ei in e) && continue
    is_QFV(Q, d, e) && return e
  end
  return "try larger"
end

