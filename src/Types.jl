########################################################################################
# Definitions of types and primitive constructors for quivers and moduli spaces
########################################################################################

"""
# Summary

`struct Quiver{T}`

A quiver is represented by its adjacency
``n \\times n`` matrix ``adjacency = (a_{ij})``,\\
where ``n`` is the number of vertices
and ``a_{ij}`` is the number of arrows ``i \\to j``.

# Fields

`adjacency :: AbstractMatrix{Int}`\\
`name      :: String`

"""
struct Quiver{T}
  adjacency::AbstractMatrix{Int}
  name::String

  """
      Quiver(adjacency, name = "")

  Construct a quiver from the adjacency matrix and optionally a name.

  # Examples

  ```jldoctest
  julia> m = [0 1; 2 0];

  julia> Quiver(m, "my quiver")
  my quiver
  ```
  """
  function Quiver(adjacency::AbstractMatrix{Int}, name::String="")
    size(adjacency, 1) != size(adjacency, 2) &&
      throw(ArgumentError("adjacency matrix must be square"))
    adj = coerce_matrix(adjacency)
    T = size(adjacency, 1)
    return new{T}(adj, name)
  end

  """
      Quiver(arrows::String)

  Construct a quiver from its arrows encoded in a string.

  The string is a comma-separated list of chains `i-j-k-...`. Within a chain, a run of
  `r` hyphens between two vertices encodes `r` arrows from the first to the second, and a
  chain `i-j-k` is read left to right, giving arrows `i → j` and `j → k`.

  Vertices are indexed as follows. If every vertex token is a positive integer, the
  integers fix the vertex indices: they are sorted and, if they are not already `1, …, n`,
  relabeled to `1, …, n` (the `k`-th smallest integer becomes vertex `k`) with a warning.
  So `"1-6,2-6,3-6,4-6,5-6,7--1"` (which uses all of `1` to `7`) keeps those indices, while
  both `"2--3"` and `"1--3"` are relabeled to `"1--2"`. Otherwise (any non-integer token)
  the tokens are merely labels, numbered `1, …, n` in order of first appearance.
  QuiverTools has no notion of vertex labels beyond this index.

  # Examples

  ```jldoctest
  julia> Quiver("1--2-3")
  Quiver with adjacency matrix [0 2 0; 0 0 1; 0 0 0]

  julia> Quiver("a---b")
  Quiver with adjacency matrix [0 3; 0 0]

  julia> Quiver("1--2,1---3,2----3")
  Quiver with adjacency matrix [0 2 3; 0 0 4; 0 0 0]
  ```

  Integer tokens index the vertices directly, so the order in the string does not matter:

  ```jldoctest
  julia> Quiver("2-1")
  Quiver with adjacency matrix [0 0; 1 0]
  ```
  """
  function Quiver(arrows::String)
    arrows = replace(arrows, r"\s" => "")
    chains = split(arrows, ",")

    # distinct vertex tokens, in order of first appearance
    tokens = String[]
    for chain in chains, token in split(chain, "-")
      tok = String(token)
      isempty(tok) || tok in tokens || push!(tokens, tok)
    end

    # if every token is a positive integer, the integers fix the vertex order; they are
    # sorted and, when not already 1, …, n, relabeled to 1, …, n (with a warning).
    # otherwise the tokens are labels, numbered 1, …, n in order of first appearance.
    ints = tryparse.(Int, tokens)
    if !isempty(tokens) && all(!isnothing, ints)
      all(>(0), ints) ||
        throw(ArgumentError("integer vertex labels must be positive"))
      labels = sort(unique(ints))
      n = length(labels)
      if labels != 1:n
        @warn "Quiver: integer vertex labels $labels are not 1:$n; " *
          "relabeling to 1:$n (the kth smallest label becomes vertex k)"
      end
      rank = Dict(v => i for (i, v) in enumerate(labels))
      index = Dict(tokens[k] => rank[ints[k]] for k in eachindex(tokens))
    else
      n = length(tokens)
      index = Dict(tok => i for (i, tok) in enumerate(tokens))
    end

    A = zeros(Int, n, n)
    for chain in chains
      pieces = split(chain, "-")
      source = index[String(pieces[1])]
      number = 1
      for piece in pieces[2:end]
        if isempty(piece)
          number += 1                  # another hyphen ⇒ one more arrow
        else
          target = index[String(piece)]
          A[source, target] += number   # accumulate parallel arrows
          number = 1
          source = target               # continue the chain from here
        end
      end
    end
    return Quiver(A, "")
  end
end

==(Q1::Quiver, Q2::Quiver) = Q1.adjacency == Q2.adjacency
hash(Q::Quiver) = hash(Q.adjacency)

function show(io::IO, Q::Quiver)
  if Q.name == ""
    print(io, "Quiver with adjacency matrix ")
    print(io, Q.adjacency)
  else
    print(io, Q.name)
  end
end

"""
# Summary

`abstract type QuiverModuli`

Abstract type for a moduli space or stack of quiver representations.

# Supertype Hierarchy

`QuiverModuliSpace <: QuiverModuli <: Any`\\
`QuiverModuliStack <: QuiverModuli <: Any`

"""
abstract type QuiverModuli end

"""
# Summary

`struct ChowRing`

A Type used to encode various Chow ring data.

# Fields

`parent :: QuiverModuli`\\
`ring   :: Singular.PolyRing{Singular.n_Q}`\\
`chi    :: AbstractVector{Int}`\\
`point  :: Union{Singular.spoly{Singular.n_Q},UndefInitializer}`\\
`todd   :: Union{Singular.spoly{Singular.n_Q},UndefInitializer}`\\
`_R     :: Singular.PolyRing{Singular.n_Q}`\\
`_inclusion :: Singular.SAlgHom{Singular.Rationals}`
"""
mutable struct ChowRing
  parent::QuiverModuli
  ring::Singular.PolyRing{Singular.n_Q}
  chi::AbstractVector{Int}
  point::Union{Singular.spoly{Singular.n_Q},UndefInitializer}
  todd::Union{Singular.spoly{Singular.n_Q},UndefInitializer}
  # `nothing` until computed; caches dimension() and lets unsafe pre-seed 1 - <d,d> (#20).
  # Only finite dimensions are cached; the empty-moduli case (`-Inf`) is left uncached, so
  # the value is a plain `Int`.
  _dimension::Union{Int,Nothing}
  _R::Singular.PolyRing{Singular.n_Q}
  _inclusion::Singular.SAlgHom{Singular.Rationals}

  """   ChowRing()"""
  function ChowRing()
    # `new()` leaves the reference fields genuinely undefined (checked via `isdefined`),
    # but a bits-union field reads as defined-with-garbage, so seed the sentinel explicitly.
    chow = new()
    chow._dimension = nothing
    return chow
  end
end

function show(io::IO, chow::ChowRing)
  ring = isdefined(chow, :ring) ? chow.ring : UndefInitializer()
  chi = isdefined(chow, :chi) ? chow.chi : UndefInitializer()
  point = isdefined(chow, :point) ? chow.point : UndefInitializer()
  todd = isdefined(chow, :todd) ? chow.todd : UndefInitializer()
  print(
    io,
    "Chow ring on
  $(chow.parent)

  Intersection theory data:

 - Chow ring: $(ring)
 - Linearization: $(chi)
 - Point class: $(point)
 - Todd class: $(todd)
    ",
  )
end

linearization(CH::ChowRing) = CH.chi

"""
# Summary

`struct QuiverModuliSpace`

The moduli space of representations of a quiver `Q` of dimension vector `d`
depends on a choice of stability parameter `theta` and on whether we consider
stable or semistable representations.

# Fields

`Q         :: Quiver`\\
`d         :: AbstractVector{Int64}`\\
`theta     :: AbstractVector{Int64}`\\
`condition :: String`\\
`denom     :: Function`

# Supertype Hierarchy

`QuiverModuliSpace <: QuiverModuli <: Any`

"""
struct QuiverModuliSpace <: QuiverModuli
  Q::Quiver
  d::AbstractVector{Int}
  theta::AbstractVector{Int}
  condition::String
  denom::Function
  chow::ChowRing
end
function QuiverModuliSpace(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  condition::String="semistable",
  denom::Function=sum,
)
  condition in ["stable", "semistable"] ||
    throw(ArgumentError("condition must be \"stable\" or \"semistable\""))
  length(d) == n_vertices(Q) ||
    throw(ArgumentError("length of d must equal the number of vertices"))
  length(theta) == n_vertices(Q) ||
    throw(ArgumentError("length of theta must equal the number of vertices"))
  d = coerce_vector(d)
  theta = coerce_vector(theta)
  M = QuiverModuliSpace(Q, d, theta, condition, denom, ChowRing())
  setfield!(M.chow, :parent, M)
  return M
end

function set_linearization!(M::QuiverModuliSpace, chi::AbstractVector{Int})
  chi' * M.d != 1 && throw(DomainError("Invalid linearization"))
  setfield!(M.chow, :chi, chi)
  return nothing
end

linearization(M::QuiverModuliSpace) = linearization(M.chow)

"""
# Summary

`struct QuiverModuliStack`

The moduli stack of representations of a quiver `Q` of dimension vector `d`
depends on a choice of stability parameter `theta` and on whether we consider
stable or semistable representations.


# Fields

`Q         :: Quiver`\\
`d         :: AbstractVector{Int64}`\\
`theta     :: AbstractVector{Int64}`\\
`condition :: String`\\
`denom     :: Function`

# Supertype Hierarchy

`QuiverModuliStack <: QuiverModuli <: Any`

"""
struct QuiverModuliStack <: QuiverModuli
  Q::Quiver
  d::AbstractVector{Int}
  theta::AbstractVector{Int}
  condition::String
  denom::Function

  function QuiverModuliStack(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int}=canonical_stability(Q, d),
    condition::String="semistable",
    denom::Function=sum,
  )
    condition in ["stable", "semistable"] ||
      throw(ArgumentError("condition must be \"stable\" or \"semistable\""))
    length(d) == n_vertices(Q) ||
      throw(ArgumentError("length of d must equal the number of vertices"))
    length(theta) == n_vertices(Q) ||
      throw(ArgumentError("length of theta must equal the number of vertices"))
    d = coerce_vector(d)
    theta = coerce_vector(theta)
    return new(Q, d, theta, condition, denom)
  end
end

function show(io::IO, M::QuiverModuliSpace)
  print(
    io,
    "Quiver moduli space defined as follows:
 - quiver: $(M.Q)
 - dimension vector: $(M.d)
 - stability parameter: $(M.theta)
 - condition: $(M.condition)
    ",
  )
end
function show(io::IO, M::QuiverModuliStack)
  print(
    io,
    "Quiver moduli stack defined as follows:
 - quiver: $(M.Q)
 - dimension vector: $(M.d)
 - stability parameter $(M.theta)
 - condition: $(M.condition)
    ",
  )
end

"""
# Summary

`struct FramedQuiverModuliSpace`

The framed quiver moduli space ``M^{\\Theta\\text{-fr}}(Q, \\mathbf{d}, \\mathbf{n})`` of a
quiver `Q`, base dimension vector `d`, base stability parameter `theta` and framing datum
`n`. It is realized as the quiver moduli space of the framed quiver ``\\widehat{Q}`` (see
[`framed_quiver`](@ref) / [`coframed_quiver`](@ref)) for dimension vector
``\\widehat{\\mathbf{d}}`` and an induced stability parameter ``\\widehat{\\Theta}``; use
[`total_space`](@ref) to obtain that ordinary [`QuiverModuliSpace`](@ref) and [`base`](@ref)
to obtain the codomain ``M^{\\Theta}(Q, \\mathbf{d})`` of the projection
``p\\colon M^{\\Theta\\text{-fr}}(Q, \\mathbf{d}, \\mathbf{n}) \\to M^{\\Theta}(Q, \\mathbf{d})``.

The fibres of `p` over a Luna stratum are described by [`fibre`](@ref); see
[arXiv:2607.12895](https://arxiv.org/abs/2607.12895).

# Fields

`Q         :: Quiver`                base quiver.\\
`d         :: AbstractVector{Int}`   base dimension vector.\\
`theta     :: AbstractVector{Int}`   base stability parameter.\\
`denom     :: Function`              slope denominator.\\
`n         :: AbstractVector{Int}`   framing datum.\\
`coframed  :: Bool`                  `false`: arrows ``i_0 \\to i`` (framing vertex first);
                                     `true`: arrows ``i \\to i_0`` (coframing vertex last).
"""
struct FramedQuiverModuliSpace
  Q::Quiver
  d::AbstractVector{Int}
  theta::AbstractVector{Int}
  denom::Function
  n::AbstractVector{Int}
  coframed::Bool
end
function FramedQuiverModuliSpace(
  Q::Quiver,
  d::AbstractVector{Int};
  n::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
  coframed::Bool=false,
)
  length(d) == n_vertices(Q) ||
    throw(ArgumentError("length of d must equal the number of vertices"))
  length(theta) == n_vertices(Q) ||
    throw(ArgumentError("length of theta must equal the number of vertices"))
  length(n) == n_vertices(Q) ||
    throw(ArgumentError("length of n must equal the number of vertices"))
  return FramedQuiverModuliSpace(
    Q, coerce_vector(d), coerce_vector(theta), denom, coerce_vector(n), coframed
  )
end

function show(io::IO, X::FramedQuiverModuliSpace)
  print(
    io,
    "$(X.coframed ? "Coframed" : "Framed") quiver moduli space defined as follows:
 - base quiver: $(X.Q)
 - base dimension vector: $(X.d)
 - base stability parameter: $(X.theta)
 - framing datum: $(X.n)
    ",
  )
end

"""
# Summary

`struct NilpotentLocus`

A marker for the closed sublocus of nilpotent representations inside a quiver moduli space
`ambient` (here a [`FramedQuiverModuliSpace`](@ref)): the representations whose underlying
representation of the quiver (ignoring the framing) is nilpotent.

The nilpotency is recorded, not computed: there is no general routine here for the
nullcone of a non-acyclic quiver, so `ambient(N)` is the ambient moduli space and the
nilpotent locus is identified by hand, as in
[arXiv:2607.12895](https://arxiv.org/abs/2607.12895). It arises as the fibre of the framed
projection, see [`fibre`](@ref).

# Fields

`ambient :: FramedQuiverModuliSpace`
"""
struct NilpotentLocus{M}
  ambient::M
end

function show(io::IO, N::NilpotentLocus)
  print(io, "Nilpotent locus inside\n  ", N.ambient)
end

"""
# Summary

`struct HNType`

A Harder-Narasimhan type for a quiver, a dimension vector `d`
and a stability parameter `theta`: an ordered tuple ``(d^1, \\dots, d^s)``
of dimension vectors with ``d^1 + \\dots + d^s = d``,
recording the dimension vectors of the semistable subquotients
of the Harder-Narasimhan filtration of a representation
of dimension vector ``d``.

The entries are ordered by strictly decreasing slope,
``\\mu(d^1) > \\dots > \\mu(d^s)``,
so that ``d^1`` is the dimension vector of the first step of the filtration,
the maximal destabilizing subrepresentation.

# Fields

 `hn :: Vector{Vector{Int}}`: the dimension vectors ``d^1, \\dots, d^s``,
 in order of strictly decreasing slope. The constructor accepts any
 `Vector{<:AbstractVector{Int}}` and coerces the entries.
"""
struct HNType
  hn::Vector{Vector{Int}}
  function HNType(dstar::Vector{<:AbstractVector{Int}})
    return new(coerce_vector.(dstar))
  end
end

function show(io::IO, H::HNType)
  print(io, "$(Vector.(H.hn))") # coercion back to vector is slow, but it's just for printing
end

==(H1::HNType, H2::HNType) = H1.hn == H2.hn
==(H::HNType, x::Vector{<:AbstractVector{Int}}) = H.hn == x
hash(H::HNType) = hash(H.hn)
length(H::HNType) = length(H.hn)
Base.getindex(H::HNType, i) = getindex(H.hn, i)
Base.iterate(H::HNType) = iterate(H.hn)
Base.iterate(H::HNType, i) = iterate(H.hn, i)
Base.getindex(H::HNType, i::Int) = getindex(H.hn, i)
Base.convert(::Type{<:HNType}, x::Vector{<:AbstractVector{Int}}) = HNType(x)

"""
# Summary

`mutable struct Bundle`

An abstract bundle on a quiver moduli. It is represented in practice by
its Chern character.

# Fields

`parent :: ChowRing`\\
`rank   :: Int`\\
`chern_character  :: Singular.spoly{Singular.n_Q}`\\
`chern_class :: Dict{Int,Singular.spoly{Singular.n_Q}}`\\
`teleman_weights :: Dict{<:HNTypes,Vector{Int}}`
"""
mutable struct Bundle
  parent::ChowRing
  rank::Int
  chern_character::Singular.spoly{Singular.n_Q}
  chern_class::Dict{Int,Singular.spoly{Singular.n_Q}}
  teleman_weights::Dict{<:HNType,Vector{Int}}
  Bundle() = new()
end

function show(io::IO, F::Bundle)
  print(
    io,
    "Bundle of rank $(rank(F))",
  )
end

function Bundle(
  parent::ChowRing, rank::Int, chern_classes::Vector{Singular.spoly{Singular.n_Q}}
)
  n = dimension(parent.parent)
  bundle = Bundle()
  setfield!(bundle, :parent, parent)
  setfield!(bundle, :rank, rank)
  cl = Dict{Int,Singular.spoly{Singular.n_Q}}(i => chern_classes[i + 1] for i in 0:n)
  setfield!(bundle, :chern_class, cl)
  return bundle
end

function Bundle(parent::ChowRing, rank::Int, chern_class::Singular.spoly{Singular.n_Q})
  n = dimension(parent.parent)
  bundle = Bundle()
  setfield!(bundle, :parent, parent)
  setfield!(bundle, :rank, rank)
  hom = __homogeneous_components(parent.parent, chern_class)
  cl = Dict{Int,Singular.spoly{Singular.n_Q}}(i => hom[i + 1] for i in 0:n)
  setfield!(bundle, :chern_class, cl)
  return bundle
end

function Bundle(parent::ChowRing, chern_character::Singular.spoly{Singular.n_Q})
  r = constant_coefficient(chern_character)
  denominator(r) != 1 && throw(DomainError("Incorrect Chern character."))
  bundle = Bundle()
  setfield!(bundle, :parent, parent)
  setfield!(bundle, :rank, Int(Singular.numerator(r)))
  setfield!(bundle, :chern_character, chern_character)
  return bundle
end

function Bundle(parent::ChowRing, chern_character::Int)
  CH = parent.ring
  bundle = Bundle()
  setfield!(bundle, :parent, parent)
  setfield!(bundle, :rank, chern_character)
  setfield!(bundle, :chern_character, CH(chern_character))
  return bundle
end

function Bundle(M::QuiverModuliSpace, chern_character::Int)
  bundle = Bundle()
  setfield!(bundle, :parent, M.chow)
  setfield!(bundle, :rank, chern_character)
  setfield!(bundle, :chern_character, M.chow.ring(chern_character))
  return bundle
end

function Bundle(M::QuiverModuliSpace, chern_character::Singular.spoly{Singular.n_Q})
  bundle = Bundle()
  setfield!(bundle, :parent, M.chow)
  r = constant_coefficient(chern_character)
  denominator(r) != 1 && throw(DomainError("Incorrect Chern character."))
  setfield!(bundle, :rank, Int(Singular.numerator(r)))
  setfield!(bundle, :chern_character, chern_character)
  return bundle
end

function Bundle(M::QuiverModuliSpace, rank::Int, x::Singular.spoly{Singular.n_Q})
  bundle = Bundle()
  setfield!(bundle, :parent, M.chow)
  setfield!(bundle, :rank, rank)
  hom = __homogeneous_components(M, x)
  cl = Dict{Int,Singular.spoly{Singular.n_Q}}(i => hom[i + 1] for i in 0:dimension(M))
  setfield!(bundle, :chern_class, cl)
  return bundle
end

function Bundle(M::QuiverModuliSpace, rank::Int, x::Dict{Int,Singular.spoly{Singular.n_Q}})
  bundle = Bundle()
  setfield!(bundle, :parent, M.chow)
  setfield!(bundle, :rank, rank)
  setfield!(bundle, :chern_class, x)
  return bundle
end

function Bundle(M::QuiverModuliSpace, rank::Int, x::Vector{Singular.spoly{Singular.n_Q}})
  bundle = Bundle()
  setfield!(bundle, :parent, M.chow)
  setfield!(bundle, :rank, rank)
  hom = __homogeneous_components(M, sum(x))
  cl = Dict{Int,Singular.spoly{Singular.n_Q}}(i => hom[i + 1] for i in 0:dimension(M))
  setfield!(bundle, :chern_class, x)
  return bundle
end

function Bundle(M::QuiverModuliSpace, weights::Dict{HNType,Vector{Int}})
  r = length(first(values(weights)))
  !all(length(v) == r for v in values(weights)) &&
    throw(ArgumentError("Incorrect weights."))
  bundle = Bundle()
  setfield!(bundle, :parent, M.chow)
  setfield!(bundle, :rank, r)
  setfield!(bundle, :teleman_weights, weights)
  return bundle
end

"""
# Summary

`struct LunaType`

A struct encoding a Luna type of a dimension vector `d` for a stability parameter `theta`.

# Description

A Luna type of `d` for `theta` is an unordered sequence
``(\\mathbf{d}^1, m_1), \\dots, (\\mathbf{d}^s, m_s)`` of pairs of dimension vectors
``\\mathbf{d}^k`` and positive integers ``m_k`` such that

- ``m_1 \\mathbf{d}^1 + \\dots + m_s \\mathbf{d}^s = \\mathbf{d}``,
- ``\\mu_\\theta(\\mathbf{d}^k) = \\mu_\\theta(\\mathbf{d})`` for all ``k``, and
- each ``\\mathbf{d}^k`` admits a ``\\theta``-stable representation.

Luna types (also called *semistable representation types* in the reference below, or
*polystable types* / *decomposition types* elsewhere) index the strata of the Luna
stratification of the moduli space ``M^{ss}_\\theta(Q, \\mathbf{d})``:
the open stratum is the stable locus, given by the trivial type `Dict(d => [1])`, and the
remaining (non-trivial) types stratify the properly semistable locus. Each stratum is
described étale-locally by a *local quiver* assembled from the stable summands; see
[`local_quiver_setting`](@ref) and
[Adriaenssens--Le Bruyn](https://mathscinet.ams.org/mathscinet/relay-station?mr=1972892).

# Encoding

A `LunaType` wraps a dictionary `data` whose keys are the (distinct) dimension vectors
``\\mathbf{d}^k`` and whose values are non-empty lists of positive integers
``[p_{k, 1}, \\dots, p_{k, t_k}]``. Such an entry encodes that ``\\mathbf{d}^k`` occurs
``t_k`` times in the sequence, coupled with the multiplicities
``p_{k, 1}, \\dots, p_{k, t_k}``, so that

```math
\\sum_k (p_{k, 1} + \\dots + p_{k, t_k})\\, \\mathbf{d}^k = \\mathbf{d}.
```

For example, `Dict([1, 1] => [2, 1])` represents the Luna type in which the dimension
vector `[1, 1]` appears twice: once with multiplicity `2` and once with multiplicity `1`
(so it contributes `(2 + 1) * [1, 1] = [3, 3]` to `d`).

The list ``[p_{k, 1}, \\dots, p_{k, t_k}]`` records ``t_k`` *distinct* stable summands of
dimension vector ``\\mathbf{d}^k``, so its length is bounded by the number of
non-isomorphic ``\\theta``-stable representations of ``\\mathbf{d}^k``: a ``\\mathbf{d}^k``
with a unique stable representation (e.g. a real Schur root) can appear only once, with a
length-one list. The three conditions above are therefore necessary but not sufficient for
a dictionary to be realized by an actual representation.

# Fields

 `data :: Dict{Vector{Int},Vector{Int}}`: dimension vectors mapped to their lists of
 multiplicities, as described above.

"""
struct LunaType
  data::Dict{Vector{Int},Vector{Int}}

  function LunaType(new_luna::Dict{<:AbstractVector{Int},Vector{Int}})
    # T = length(collect(keys(new_luna))[1])

    return new(Dict(coerce_vector(tau) => new_luna[tau] for tau in keys(new_luna)))
  end
end

function show(io::IO, L::LunaType)
  print(io, "Dict(")
  dim_vector = collect(keys(L.data))
  l = length(dim_vector)
  for i in 1:(l - 1)
    print(io, "$(Vector(dim_vector[i])) => $(L.data[dim_vector[i]]), ")
  end
  print(io, "$(Vector(dim_vector[l])) => $(L.data[dim_vector[l]]))")
end

==(L1::LunaType, L2::LunaType) = L1.data == L2.data
==(L::LunaType, x::Dict{<:AbstractVector{Int},Vector{Int}}) = L.data == x
hash(L::LunaType) = hash(L.data)
length(L::LunaType) = length(L.data)
getindex(L::LunaType, i) = getindex(L.data, i)
setindex!(L::LunaType, value, key...) = setindex!(L.data, value, key...)
iterate(L::LunaType) = iterate(L.data)
iterate(L::LunaType, i) = iterate(L.data, i)
keys(L::LunaType) = keys(L.data)
haskey(L::LunaType, key) = haskey(L.data, key)
