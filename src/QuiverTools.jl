module QuiverTools

using Pkg
using StaticArrays

using Memoization: Memoization
using IterTools: IterTools
using LinearAlgebraX: LinearAlgebraX
using Singular: Singular
using AbstractAlgebra: AbstractAlgebra
using Nemo: Nemo

import Base.show, Base.==, Base.hash, Base.getindex, Base.length, Base.iterate, Base.keys
import Base: *, +, -, ^
import Memoization: @memoize
import IterTools: subsets
import LinearAlgebraX: rankx
import Singular:
  polynomial_ring,
  degree,
  coeff,
  constant_coefficient,
  AlgebraHomomorphism,
  preimage,
  Ideal,
  quotient_ideal,
  QuotientRing,
  fraction_field,
  std,
  gens,
  base_ring
import Combinatorics: combinations, with_replacement_combinations, partitions

export n_vertices,
  n_arrows, arrows, indegree, outdegree, is_acyclic, is_connected, is_sink, is_source
export euler_form, canonical_stability, is_coprime, slope
export underlying_graph, euler_matrix
export is_schur_root, is_real_root, is_imaginary_root, is_isotropic_root,
  generic_ext, generic_hom, canonical_decomposition, in_fundamental_domain
export all_hn_types,
  is_hn_type, has_semistables, has_stables, codimension_hn_stratum, is_amply_stable
export is_generic_subdimension_vector, all_generic_subdimension_vectors

# TODO add missing doctests across codebase
# TODO add safety checks everywhere in the codebase

include("Types.jl")

import Pkg

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])

function _print_banner()
  printstyled(raw"""   ___"""; color=:red)
  printstyled(raw"""       _             """)
  println("  |")
  printstyled(raw"""  / _ \ """; color=:red)
  printstyled(raw"""_  _(_)_ _____ _ _ """)
  println("  |  Software package for quivers")
  printstyled(raw""" | (_) | """; color=:red)
  printstyled(raw"""|| | \ V / -_) '_|""")
  println("  |  and moduli of their representations")
  printstyled(raw"""  \__\_\\"""; color=:red)
  printstyled(raw"""\_,_|_|\_/\___|_|  """)
  println("  |")
  printstyled(raw"""       _____         _    """; color=:yellow)
  println("   |")
  printstyled(raw"""      |_   _|__  ___| |___"""; color=:yellow)
  println("   |  Manual: https://julia.quiver.tools")
  printstyled(raw"""        | |/ _ \/ _ \ (_-<"""; color=:yellow)
  println("   |  Version $(VERSION_NUMBER)")
  printstyled(raw"""        |_|\___/\___/_/__/"""; color=:yellow)
  return println("   |")
end

function __init__()
  if displaysize(stdout)[2] >= 80
    _print_banner()
  end

  return nothing
end

# TODO move all quiver things to Quiver.jl, so that QuiverTools.jl is only meta things

function deglex_key(Q::Quiver, e::AbstractVector{Int})::Int
  b = maximum(e) + 1
  n = n_vertices(Q)

  return (sum(e[i] * b^(n - i) for i in 1:length(e)) + sum(e) * b^n)
end

"""
    underlying_graph(Q::Quiver)

Returns the (necessarily symmetric) adjacency matrix
of the underlying graph of the quiver.

# Input

- `Q::Quiver` a quiver

# Output

- the adjacency matrix of the underlying graph of the quiver

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> underlying_graph(Q) == [0 4; 4 0]
true
```
"""
function underlying_graph(Q::Quiver)
  return Matrix{Int}(Q.adjacency + transpose(Q.adjacency) - diagonal(Q.adjacency))
end

"""
    n_vertices(Q::Quiver)

Returns the number of vertices of the quiver.

# Input

- `Q::Quiver` a quiver

# Output

- the number of vertices of the quiver

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> n_vertices(Q) == 2
true
```
"""
n_vertices(Q::Quiver) = size(Q.adjacency)[1]

"""
    n_arrows(Q::Quiver)

Returns the number of arrows of the quiver.

# Input

- `Q::Quiver` a quiver

# Output

- the number of arrows of the quiver

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> n_arrows(Q) == 4
true
```
"""
n_arrows(Q::Quiver) = sum(Q.adjacency)

"""
    is_acyclic(Q::Quiver)

Checks whether the quiver is acyclic, i.e. has no oriented cycles.

# Input

- `Q::Quiver` a quiver

# Output

- `true` if the quiver is acyclic, `false` otherwise

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> is_acyclic(Q)
true
```
"""
is_acyclic(Q::Quiver) = all(entry == 0 for entry in Q.adjacency^n_vertices(Q))

"""
    is_connected(Q::Quiver)

Checks whether the quiver is connected.

# Input

- `Q::Quiver` a quiver

# Output

- `true` if the underlying graph of the quiver is connected, `false` otherwise

# Examples

```jldoctest
julia> Q = Quiver([0 1 0; 0 0 1; 1 0 0]);

julia> is_connected(Q)
true

julia> Q = Quiver([0 1 0; 1 0 0; 0 0 2]);

julia> is_connected(Q)
false

julia> # The 4-Kronecker quiver:

julia> Q = kronecker_quiver(4);

julia> is_connected(Q)
true

julia> # The 4-loop quiver:

julia> Q = loop_quiver(4);

julia> is_connected(Q)
true

julia> # The 4-subspace quiver:

julia> Q = subspace_quiver(4);

julia> is_connected(Q)
true
```
"""
function is_connected(Q::Quiver)
  paths = underlying_graph(Q)
  for i in 2:(n_vertices(Q) - 1)
    paths += paths * underlying_graph(Q)
  end
  for i in 1:n_vertices(Q), j in 1:n_vertices(Q)
    if i != j && paths[i, j] == 0 && paths[j, i] == 0
      return false
    end
  end
  return true
end

"""
    indegree(Q::Quiver, j::Int)

Returns the number of incoming arrows to the vertex `j`.

# Input

- `Q::Quiver` a quiver
- `j::Int`: an integer representing the vertex

# Output

- the number of incoming arrows to the vertex `j`

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> indegree(Q, 1)
0

julia> indegree(Q, 2)
4
```
"""
indegree(Q::Quiver, j::Int) = sum(Q.adjacency[:, j])

"""
    outdegree(Q::Quiver, i::Int)

Returns the number of outgoing arrows from the vertex `i`.

# Input

- `Q::Quiver` a quiver
- `i::Int`: an integer representing the vertex

# Output

- the number of outgoing arrows from the vertex `i`

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> outdegree(Q, 1)
4

julia> outdegree(Q, 2)
0
```
"""
outdegree(Q::Quiver, i::Int) = sum(Q.adjacency[i, :])

"""
    is_source(Q::Quiver, i::Int)

Checks if the vertex `i` is a source, i.e., a vertex with no incoming arrows.

# Input

- `Q::Quiver` a quiver
- `i::Int`: an integer representing the vertex

# Output

- `true` if the vertex `i` is a source, `false` otherwise

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> is_source(Q, 1)
true

julia> is_source(Q, 2)
false
```
"""
is_source(Q::Quiver, i::Int) = indegree(Q, i) == 0

"""
    is_sink(Q::Quiver, j::Int)

Checks if the vertex `j` is a sink, i.e., a vertex with no outgoing arrows.

# Input

- `Q::Quiver` a quiver
- `j::Int`: an integer representing the vertex

# Output

- `true` if the vertex `j` is a sink, `false` otherwise

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> is_sink(Q, 1)
false

julia> is_sink(Q, 2)
true
```
"""
is_sink(Q::Quiver, j::Int) = outdegree(Q, j) == 0

"""
    arrows(Q::Quiver)

Returns a list of all arrows of the quiver `Q`.

# Input

- `Q::Quiver` a quiver

# Output

- a list of all arrows of the quiver `Q`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> arrows(Q)
3-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 2]
 [1, 2]
```
"""
function arrows(Q::Quiver)
  n = n_vertices(Q)
  return reduce(
    vcat,
    [[i, j] for k in 1:Q.adjacency[i, j]] for i in 1:n for
    j in 1:n if Q.adjacency[i, j] > 0
  )
end

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

"""
    euler_matrix(Q::Quiver)

Returns the Euler matrix of the quiver.

The Euler matrix of a quiver ``Q`` is defined as
```math
E = I - A,
```
where ``A`` is the adjacency matrix of ``Q`` and ``I``
is the identity matrix of the same size as ``A``.

# Input

- `Q::Quiver` a quiver

# Output

- the Euler matrix of the quiver

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> euler_matrix(Q) == [1 -4; 0 1]
true
```
"""
@memoize Dict euler_matrix(Q::Quiver) = identity_matrix(n_vertices(Q)) - Q.adjacency

"""
    euler_form(Q::Quiver, x, y)

Computes the Euler form of the quiver for vectors `x` and `y`.

The Euler form is defined as the bilinear form
```math
\\langle x,y\\rangle = x^T * E * y,
```
where ``E`` is the Euler matrix of the quiver.

# Input

- `Q::Quiver` a quiver
- `x::AbstractVector{Int}`: a vector
- `y::AbstractVector{Int}`: a vector

# Output

- the Euler form ``\\langle x, y\\rangle_{Q}``

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> euler_form(Q, [1, 1], [1, 1]) == -2
true
```
"""
euler_form(Q::Quiver, x::AbstractVector{Int}, y::AbstractVector{Int}) =
  x' * euler_matrix(Q) * y

"""
    canonical_stability(Q::Quiver, d::AbstractVector{Int})

The canonical stability parameter for the couple ``(Q, d)`` is given by ``<d,-> - <-,d>``

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector

# Output

- the canonical stability parameter for `Q` and `d`

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2,3];

julia> canonical_stability(Q, d) == [9, -6]
true
```
"""
function canonical_stability(Q::Quiver, d::AbstractVector{Int})
  return -(-transpose(euler_matrix(Q)) + euler_matrix(Q)) * d
end

"""
    is_coprime(d, theta)

Checks whether the given dimension vector `d` is `theta`-coprime for
the stability parameter ``\\theta``.

A dimension vector ``d`` is said to be ``\\theta``-coprime for the
stability parameter ``\\theta`` if all subdimension vectors ``0 \\neq e < d``
satisfy ``\\theta * e \\neq 0``.

# Input

- `d::AbstractVector{Int}` a dimension vector
- `theta::AbstractVector{Int}` a stability parameter

# Output

- `true` if `d` is `theta`-coprime, `false` otherwise

# Examples

```jldoctest
julia> d = [2, 3]; theta = [3, -2];

julia> is_coprime(d, theta)
true

julia> is_coprime([3, 3], theta)
false
```
"""
function is_coprime(d::AbstractVector{Int}, theta::AbstractVector{Int})
  return all(
    e -> theta' * e != 0,
    all_subdimension_vectors(d; nonzero=true, strict=true),
  )
end

"""
    is_coprime(d)

Checks if the gcd of all the entries of d is ``1``.

# Input

- `d::AbstractVector{Int}` a vector

# Output

- `true` if the gcd of all the entries of `d` is ``1``, `false` otherwise

# Examples

```jldoctest
julia> is_coprime([2, 3])
true

julia> is_coprime([3, 3])
false
```
"""
is_coprime(d::AbstractVector{Int}) = gcd(d) == 1

"""
    slope(d, theta, denom=sum)

Returns the slope of `d`
with respect to the stability parameter ``\\theta``
and a choice of a denominator function `denom`.

The slope function for ``\\theta`` and ``\\alpha`` is defined as
```math
\\mu = \\frac{\\theta}{\\alpha} := x \\mapsto \\frac{\\theta \\cdot x}{\\alpha(x)}.
```

# Input

- `d::AbstractVector{Int}` a dimension vector
- `theta::AbstractVector{Int}` a stability parameter
- `denom::Function` a function to compute the denominator. Default is `sum`

# Output

- the slope of `d` with respect to `theta` and `denom`

# Examples

```jldoctest
julia> slope([2,3], [3,-2])
0//1
```
"""
function slope(d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
  return (theta' * d)//denom(d)
end

"""
    all_destabilizing_subdimension_vectors(d, theta, denom=sum)

Returns the subdimension vectors of `d` with a strictly larger slope than `d`.

# Input

- `d::AbstractVector{Int}` a dimension vector
- `theta::AbstractVector{Int}` a stability parameter
- `denom::Function` a function to compute the denominator. Default is `sum`

# Output

- an array of subdimension vectors of `d` with a strictly larger slope than `d`
"""
@memoize Dict function all_destabilizing_subdimension_vectors(
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  # as silly as it looks this is faster.
  # TODO faster than what? is the calculation `b = slope(d, theta, denom)` not being reused in the lambda expression?
  # if that's the point, maybe it makes sense to cache the lambda function instead of the result of `slope(d, theta, denom)`?
  b = slope(d, theta, denom)
  return filter(
    e -> slope(e, theta, denom) > b,
    all_subdimension_vectors(d; nonzero=true),
  )
end

"""
    has_semistables(Q::Quiver, d, theta=canonical_stability(Q, d), denom=sum)

Checks if `Q` admits a ``\\theta``-semistable representation of dimension vector ``d``.

A representation ``V`` is said to be ``\\mu``-semistable if
for all of its subrepresentations ``W`` with ``dim(W) < dim(V)``, we have
```math
\\mu(\\dim(W)) \\leq \\mu\\dim((V)).
```

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector
- `theta::AbstractVector{Int}` a stability parameter. Default is `canonical_stability(Q, d)`
- `denom::Function` a function to compute the denominator. Default is `sum`

# Output

- `true` if there is a `theta`-semistable representation of dimension vector `d`,
  `false` otherwise

# Examples

```jldoctest
julia> A2 = kronecker_quiver(1); theta = [1,-1];

julia> has_semistables(A2, [1,1], theta)
true

julia> has_semistables(A2, [2,2], theta)
true

julia> has_semistables(A2, [1,2], theta)
false

julia> has_semistables(A2, [0,0], theta)
true
```
The 3-Kronecker quiver:

```jldoctest
julia> K3 = kronecker_quiver(3); theta = [3,-2];

julia> has_semistables(K3, [2,3], theta)
true

julia> has_semistables(K3, [1,4], theta)
false
```
"""
@memoize Dict function has_semistables(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
)
  # TODO Julia shorthand? + explain this is a base case
  if all(di == 0 for di in d)
    return true
  else
    # collect the list of all subdimension vectors e of bigger slope than d
    slope_d = slope(d, theta, denom)
    # TODO this variable doesn't follow our conventions, why not `subdimensions_bigger_slope`?
    subdimensionsBiggerSlope = filter(
      e -> slope(e, theta, denom) > slope_d,
      all_subdimension_vectors(d; nonzero=true, strict=true),
    )
    # to have semistable representations, none of the vectors above must be
    # a generic subdimension vector.
    return all(e -> !is_generic_subdimension_vector(Q, e, d), subdimensionsBiggerSlope)
  end
end

"""
    has_stables(Q::Quiver, d, theta=canonical_stability(Q, d), denom=sum)

Checks if Q admits a ``theta``-stable representation of dimension vector ``d``.

A representation ``V`` is said to be ``mu``-stable if
for all of its subrepresentations ``W`` with ``dim(W) < dim(V)``, we have
```math
\\mu(\\dim(W)) < \\mu\\dim((V)).
```

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector
- `theta::AbstractVector{Int}` a stability parameter. Default is `canonical_stability(Q, d)`
- `denom::Function` a function to compute the denominator. Default is `sum`

# Output

- `true` if there is a `theta`-stable representation of dimension vector `d`,
  `false` otherwise

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2, 3]; theta = [3, -2];

julia> has_stables(Q, d, theta)
true

julia> Q = kronecker_quiver(2); d = [2,2]; theta = [1,-1];

julia> has_stables(Q, d, theta)
false

julia> has_semistables(Q, d, theta)
true
```

The zero dimension vector has no stables:
```jldoctest
julia> Q = kronecker_quiver(3); d = [0,0];

julia> has_stables(Q, d)
false
```
"""
@memoize Dict function has_stables(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
)
  # TODO Julia shorthand?
  if all(di == 0 for di in d)
    return false
  else
    # collect the list of all subdimension vectors e of bigger slope than d
    slope_d = slope(d, theta, denom)
    subdimensions_bigger_or_equal_slope = filter(
      e -> slope(e, theta, denom) >= slope_d,
      all_subdimension_vectors(d; nonzero=true, strict=true),
    )
    # to have semistable representations,
    # none of the vectors above must be generic subdimension vectors.
    return all(
      e -> !is_generic_subdimension_vector(Q, e, d),
      subdimensions_bigger_or_equal_slope,
    )
  end
end

# TODO the cited paper is published
"""
    is_schur_root(Q::Quiver, d)

Checks if ``d`` is a Schur root for ``Q``.

By [Lemma 4.2, arXiv:0802.2147](https://doi.org/10.48550/arXiv.0802.2147),
this is equivalent to the existence of a stable representation of dimension vector ``d``
for the canonical stability parameter.

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector

# Output

- `true` if `d` is a Schur root for `Q`, `false` otherwise

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2,3];

julia> is_schur_root(Q, d)
true
```
"""
is_schur_root(Q::Quiver, d::AbstractVector{Int}) =
  has_stables(Q, d, canonical_stability(Q, d))

"""
    is_real_root(Q::Quiver, d)

Checks whether `d` is a real root, i.e., if ``<d, d> = 1``.
"""
is_real_root(Q, d) = euler_form(Q, d, d) == 1

"""
    is_imaginary_root(Q::Quiver, d)

Checks whether `d` is an imaginary root, i.e., if ``<d, d> \\geq 0``.
"""
is_imaginary_root(Q, d) = euler_form(Q, d, d) <= 0

"""
    is_isotropic_root(Q::Quiver, d)

Checks whether `d` is an isotropic root, i.e., if ``<d, d> = 0``.
"""
is_isotropic_root(Q, d) = euler_form(Q, d, d) == 0

"""
    is_generic_subdimension_vector(Q::Quiver, e, d)

Checks if ``e`` is a generic subdimension vector of ``d``.

A dimension vector ``e`` is called a generic subdimension vector of ``d``
if a generic representation of dimension vector ``d`` possesses a subrepresentation
of dimension vector ``e``.

By [Theorem 5.3, arXiv:0802.2147](https://doi.org/10.48550/arXiv.0802.2147),
``e`` is a generic subdimension vector of ``d`` if and only if
```math
<e',d-e> \\geq 0
```
for all generic subdimension vectors ``e'`` of ``e``.

# Input

- `Q::Quiver` a quiver
- `e::AbstractVector{Int}` a dimension vector
- `d::AbstractVector{Int}` a dimension vector

# Output

- `true` if `e` is a generic subdimension vector of `d`, `false` otherwise

# Examples

Trivial examples on the 3-Kronecker quiver:

```jldoctest
julia> Q = kronecker_quiver(3); e = [1, 2]; d = [2, 3];

julia> is_generic_subdimension_vector(Q, e, d)
true

julia> is_generic_subdimension_vector(Q, [0, 0], d)
true

julia> is_generic_subdimension_vector(Q, [2, 3], d)
true

julia> is_generic_subdimension_vector(Q, [2, 1], d)
false
```
"""
@memoize Dict function is_generic_subdimension_vector(
  Q::Quiver,
  e::AbstractVector{Int},
  d::AbstractVector{Int},
)::Bool # TODO is specifying the output type a good idea? https://docs.julialang.org/en/v1/manual/functions/#The-return-Keyword says they are not a good idea usually
  # TODO Julia shorthand?
  if e == d || all(ei == 0 for ei in e)
    return true
  end
  # to speed up computation of <eprime,d-e>
  partial_evaluation = euler_matrix(Q) * (d - e)
  # considering subdimension vectors that violate the numerical condition
  subdimensions = filter(
    eprime -> eprime' * partial_evaluation < 0, all_subdimension_vectors(e)
  )
  # none of the subdimension vectors violating the condition should be generic
  return all(eprime -> !is_generic_subdimension_vector(Q, eprime, e), subdimensions)
  # return generic_ext(Q, e, d - e) == 0 # TODO test performance
end

"""
    all_generic_subdimension_vectors(Q::Quiver, d)

Returns the list of all generic subdimension vectors of ``d``.

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector

# Output

- a list of all generic subdimension vectors of `d`

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> QuiverTools.all_generic_subdimension_vectors(Q, [2, 3])
7-element Vector{StaticArraysCore.SVector{2, Int64}}:
 [0, 0]
 [0, 1]
 [0, 2]
 [1, 2]
 [0, 3]
 [1, 3]
 [2, 3]

julia> QuiverTools.all_generic_subdimension_vectors(Q, [3, 0])
4-element Vector{StaticArraysCore.SVector{2, Int64}}:
 [0, 0]
 [1, 0]
 [2, 0]
 [3, 0]
```
"""
@memoize Dict function all_generic_subdimension_vectors(Q::Quiver, d::AbstractVector{Int})
  return filter(e -> is_generic_subdimension_vector(Q, e, d), all_subdimension_vectors(d))
end

"""
    all_hn_types(Q::Quiver, d, theta, denom=sum; unstable::Bool=false, ordered::Bool=false)

Returns a list of all the Harder--Narasimhan types of representations of ``Q``
with dimension vector ``d``, with respect to the slope function theta/denom.

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector
- `theta::AbstractVector{Int}` a stability parameter
- `denom::Function` a function to compute the denominator. Default is `sum`

Keyword inputs:

- `unstable`: if `true` exclude the trivial Harder--Narasimhan type (d),
which corresponds to stable representations. Default is `false`
- `ordered`: if `true` return the list of all HN types in ascending order.
Default is `true`

# Output

- a list of all the Harder--Narasimhan types
of representations of `Q` with dimension vector `d`

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> all_hn_types(Q, d, theta; ordered=true)
8-element Vector{HNType}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> all_hn_types(Q, [3,0], [0,0]) == [[[3, 0]]]
true

julia> Q = three_vertex_quiver(1, 4, 1); d = [4, 1, 4];

julia> theta = canonical_stability(Q, d);

julia> all_hn_types(Q, d, theta; ordered=true)
106-element Vector{HNType}:
 [[4, 1, 4]]
 [[4, 1, 3], [0, 0, 1]]
 [[4, 0, 3], [0, 1, 1]]
 [[4, 0, 3], [0, 1, 0], [0, 0, 1]]
 [[3, 1, 2], [1, 0, 2]]
 [[3, 1, 2], [1, 0, 1], [0, 0, 1]]
 [[3, 0, 2], [1, 1, 2]]
 [[3, 0, 2], [0, 1, 0], [1, 0, 2]]
 [[3, 0, 2], [1, 0, 1], [0, 1, 1]]
 [[3, 0, 2], [1, 1, 1], [0, 0, 1]]
 ⋮
 [[3, 0, 0], [1, 1, 2], [0, 0, 2]]
 [[3, 0, 0], [0, 1, 0], [1, 0, 4]]
 [[3, 0, 0], [0, 1, 0], [1, 0, 3], [0, 0, 1]]
 [[3, 0, 0], [0, 1, 0], [1, 0, 2], [0, 0, 2]]
 [[3, 0, 0], [1, 0, 1], [0, 1, 1], [0, 0, 2]]
 [[3, 0, 0], [1, 1, 1], [0, 0, 3]]
 [[3, 0, 0], [1, 1, 0], [0, 0, 4]]
 [[4, 0, 0], [0, 1, 1], [0, 0, 3]]
 [[4, 0, 0], [0, 1, 0], [0, 0, 4]]
```
"""
@memoize Dict function all_hn_types(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum;
  unstable::Bool=false,
  ordered::Bool=true,
)
  # TODO Julia shorthand?
  if all(di == 0 for di in d)
    return [HNType([zero_vector(Q)])]
  end
  # We consider just proper subdimension vectors which admit a semistable
  # representation and for which μ(e) > μ(d)
  # Note that we also eliminate d by the following
  subdimensions = filter(
    e -> has_semistables(Q, e, theta, denom),
    all_destabilizing_subdimension_vectors(d, theta, denom),
  )

  # We sort the subdimension vectors by slope because that will return the list of
  # all HN types in ascending order with respect to the partial order from
  # Definition 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
  if ordered
    subdimensions = sort(subdimensions; by=e -> slope(e, theta, denom))
  end

  # The HN types which are not of the form (d) are (e,f^1,...,f^s) where e is a
  # proper semistable subdimension vector with μ(e) > μ(d), (f^1,...,f^s) is a HN
  # type of f = d-e and μ(e) > μ(f^1) holds.

  alltypes = HNType[
    HNType(vcat([e], efstar.hn))

    for e in subdimensions for efstar in filter(
      fstar -> slope(e, theta, denom) > slope(fstar[1], theta, denom),
      all_hn_types(Q, d - e, theta, denom; ordered=ordered),
    )
  ]

  # Possibly add d again, at the beginning, because it is smallest
  # with respect to the partial order from Definition 3.6
  if !unstable && has_semistables(Q, d, theta, denom)
    pushfirst!(alltypes, HNType([d]))
  end
  return alltypes
end

"""
	is_hn_type(Q::Quiver, d, dstar, theta, denom=sum)

Checks if the given ordered list of subdimension vectors ``dstar`` is an HN type
for the datum ``(Q, d)`` and the slope stability given by ``(theta, denom)``.

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector
- `dstar::HNType` an Harder--Narasimhan type
- `theta::AbstractVector{Int}` a stability parameter
- `denom::Function` a function to compute the denominator. Default is `sum`

# Output

- `true` if `dstar` is an HN type for `Q`, `d` and the slope `theta`/`denom`, `false` otherwise

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> d = [2, 3]; dstar = [d];

julia> is_hn_type(Q, d, dstar)
true
```
"""
function is_hn_type(
  Q::Quiver,
  d::AbstractVector{Int},
  dstar::HNType,
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
)::Bool
  sum(dstar) != d && throw(ArgumentError("$(dstar) does not sum to $(d)."))

  if !all(
    slope(dstar[i], theta, denom) > slope(dstar[i + 1], theta, denom) for
    i in 1:(length(dstar) - 1)
  )
    return false
  end

  if !all(has_semistables(Q, dstari, theta, denom) for dstari in dstar)
    return false
  end
  return true
end
is_hn_type(Q::Quiver,
  d::AbstractVector{Int},
  dstar::Vector{<:AbstractVector{Int}};
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
)::Bool =
  is_hn_type(Q, d, HNType(dstar), theta, denom)

"""
    codimension_hn_stratum(Q::Quiver, stratum)

Returns the codimension of the given HN stratum.

# Input

- `Q::Quiver` a quiver
- `stratum::HNType`: an HN type

# Output

- the codimension of the HN stratum as an integer

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> HN = all_hn_types(Q, d, theta; ordered=true);

julia> [codimension_hn_stratum(Q, stratum) for stratum in HN]
8-element Vector{Int64}:
  0
  3
  4
 10
  8
  9
 12
 18
```
"""
function codimension_hn_stratum(Q::Quiver, stratum::HNType)
  length(stratum) == 1 && return 0

  return -sum(
    euler_form(Q, stratum[i], stratum[j])
    for i in 1:(length(stratum) - 1)
    for j in (i + 1):length(stratum)
  )
end
codimension_hn_stratum(Q::Quiver, stratum::Vector{<:AbstractVector{Int}}) =
  codimension_hn_stratum(Q, HNType(stratum))

"""
    is_amply_stable(Q::Quiver, d, theta, denom=sum)

Checks whether the dimension vector ``d`` is amply stable
with respect to the slope function `theta`/`denominator`.

This means that the codimension of the unstable locus
in the parameter space is at least ``2``.

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector
- `theta::AbstractVector{Int}` a stability parameter
- `denom::Function` a function to compute the denominator. Default is `sum`

# Output

- `true` if `d` is amply stable, `false` otherwise

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); d = [2, 3];

julia> is_amply_stable(Q, d, [3, -2])
true

julia> is_amply_stable(Q, d, [-3, 2])
false

julia> is_amply_stable(Q, [3, 0], [0, -3])
true
```
"""
function is_amply_stable(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  hn_types = all_hn_types(Q, d, theta, denom; unstable=true)
  return all(stratum -> codimension_hn_stratum(Q, stratum) >= 2, hn_types)
end

########################################################################################
# Canonical decomposition
########################################################################################

"""
    generic_ext(Q::Quiver, a, b)

Computes the dimension of the ``\\mathrm{Ext}^1`` group between generic representations
of dimension vectors ``a`` and ``b``.

According to [Theorem 5.4, MR1162487]
(https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487),
we have

```math
ext(a,b)=max\\{-\\langle c,b\\rangle~~|~~c~\\text{is a generic subdimension vector of }a\\}.
```

# Input

- `Q::Quiver` a quiver
- `a::AbstractVector{Int}`: a vector
- `b::AbstractVector{Int}`: a vector

# Output

- the dimension of the generic extensions ``\\mathrm{ext}^1(a, b)``

# Examples

```jldoctest
julia> Q1 = kronecker_quiver(3);

julia> generic_ext(Q1, [2, 3], [6, 7])
9

julia> generic_ext(Q1, [1, 1], [1, 0])
0

julia> Q2 = three_vertex_quiver(1, 6, 7);

julia> generic_ext(Q2, [5, 6, 7], [6, 7, 8])
483
```
"""
function generic_ext(Q::Quiver, a::AbstractVector{Int}, b::AbstractVector{Int})
  return maximum(-euler_form(Q, c, b) for c in all_generic_subdimension_vectors(Q, a))
end

"""
    generic_hom(Q::Quiver, a, b)

Computes the dimension of the ``\\mathrm{Hom}`` group between generic representations
of dimension vectors ``a`` and ``b``.

# Input

- `Q::Quiver` a quiver
- `a::AbstractVector{Int}`: a vector
- `b::AbstractVector{Int}`: a vector

# Output

- the dimension of the generic homomorphisms ``\\mathrm{hom}(a, b)``

# Examples

```jldoctest
julia> Q1 = kronecker_quiver(3);

julia> generic_hom(Q1, [2, 3], [6, 7])
0

julia> generic_hom(Q1, [1, 1], [1, 0])
1

julia> Q2 = three_vertex_quiver(1, 6, 7);

julia> generic_hom(Q2, [5, 6, 7], [6, 7, 8])
0
```
"""
function generic_hom(Q::Quiver, a::AbstractVector{Int}, b::AbstractVector{Int})
  return euler_form(Q, a, b) + generic_ext(Q, a, b)
end

"""
    canonical_decomposition(Q::Quiver, d)

Computes the canonical decomposition of the dimension vector ``d``
for the given quiver ``Q``.

If ``\\beta_1, \\dots, \\beta_{\\ell}`` is a sequence
of Schur roots such that, for all ``i \\neq j``, one has

```math
\\mathrm{ext}(\\beta_i, \\beta_j) = \\mathrm{ext}(\\beta_j, \\beta_i) = 0,
```

then the general representation of dimension ``\\sum_i \\beta_i`` is
isomorphic to the direct sum of irreducible representations
of dimension vectors ``\\beta_i``.

Such a decomposition is called the canonical decomposition.

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector

# Output

- a list of dimension vectors representing the canonical decomposition of `d`

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> canonical_decomposition(Q, [6, 7]) == [[6, 7]]
true

julia> canonical_decomposition(Q, [1, 1]) == [[1, 1]]
true

julia> canonical_decomposition(Q, [6, 2]) == [[3, 1], [3, 1]]
true

julia> Q = kronecker_quiver(2);

julia> canonical_decomposition(Q, [8, 8]) == [[1, 1] for i in 1:8]
true
```
"""
function canonical_decomposition(Q::Quiver, d::AbstractVector{Int})
  generic_subdimension_vectors = filter(e -> e != d, all_generic_subdimension_vectors(Q, d))
  for e in generic_subdimension_vectors
    if d - e in generic_subdimension_vectors &&
      generic_ext(Q, e, d - e) == 0 &&
      generic_ext(Q, d - e, e) == 0
      return vcat(canonical_decomposition(Q, e), canonical_decomposition(Q, d - e))
    end
  end
  return [d] # if nothing above worked then d is a Schur root.
end

"""
    in_fundamental_domain(Q::Quiver, d; interior::Bool=false)

Checks if the dimension vector ``d`` is in the fundamental domain of the quiver ``Q``.

The fundamental domain is the cone of dimension vectors in ``\\mathbb{Z}^{Q_0}``
such that the symmetric Tits form is negative on all the simple roots, i.e.,
for all vertices i,

```math
(s_i, d) := \\langle d, s_i\\rangle + \\langle s_i, d\\rangle  \\leq 0,
```

where ``s_i`` is the dimension vector with all entries set to ``0`` and the i-th
set to ``1``.

# Input

- `Q::Quiver` a quiver
- `d::AbstractVector{Int}` a dimension vector

Keyword arguments:

- `interior`: if `true` checks whether `d` belongs to the interior of the fundamental
domain. Default is `false`

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> in_fundamental_domain(Q, [2, 3])
true

julia> in_fundamental_domain(Q, [1, 1])
true

julia> in_fundamental_domain(Q, [2, 2])
true

julia> in_fundamental_domain(Q, [1, 2])
false
```
"""
function in_fundamental_domain(Q::Quiver, d::AbstractVector{Int}; interior::Bool=false)
  # https://arxiv.org/abs/2209.14791 uses a strict inequality,
  # while https://arxiv.org/abs/2310.15927 uses a non-strict.
  # here we set it to non-strict by default.

  simples = [unit_vector(n_vertices(Q), i) for i in 1:n_vertices(Q)]
  if interior
    return all(
      simple -> euler_form(Q, d, simple) + euler_form(Q, simple, d) < 0,
      simples,
    )
  end
  return all(simple -> euler_form(Q, d, simple) + euler_form(Q, simple, d) <= 0, simples)
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

#######################################################
# Include all the submodules
#######################################################

include("Constructors.jl")
include("Moduli.jl")
include("Teleman.jl")
include("Bundles.jl")

######################
# end of QuiverTools
######################
end
