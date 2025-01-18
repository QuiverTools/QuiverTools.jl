########################################################################################
# Definitions of types and primitive constructors for quivers and moduli spaces
########################################################################################
import Base.getindex, Base.length, Base.iterate, Base.hash, Base.==, Base.show
export Quiver, HNType, QuiverModuli, QuiverModuliSpace, QuiverModuliStack, Bundle

"""
# Summary

`struct Quiver`

A quiver is represented by its adjacency
``n \\times n`` matrix ``adjacency = (a_{ij})``,\\
where ``n`` is the number of vertices
and ``a_{ij}`` is the number of arrows ``i \\to j``.

# Fields

`adjacency :: AbstractMatrix{Int}`\\
`name      :: String`

"""
struct Quiver
  adjacency::AbstractMatrix{Int}
  name::String

  """
      Quiver(adjacency, name = "")

  Constructs a quiver starting from its adjacency matrix, and an optional name.

  # Examples

  ```jldoctest
  julia> m = [0 1; 2 0];

  julia> Quiver(m, "my quiver")
  my quiver, with adjacency matrix [0 1; 2 0]
  ```
  """
  function Quiver(adjacency::AbstractMatrix{Int}, name::String="")
    if !(size(adjacency)[1] == size(adjacency)[2])
      throw(DomainError(adjacency, "adjacency matrix must be square"))
    else
      adj = SMatrix{size(adjacency)...}(adjacency)
      new(adj, name)
    end
  end

  """
      Quiver(arrows)

  Constructs a quiver based on its arrows encoded in a string.

  the string `arrows` must be of the form

  ```i---j,k-...-s```

  where `i`, `j` and all the vertices are positive integers.
  The amount of characters between `i` and  `j` is then the number of arrows i -> j.

  # Examples

  ```jldoctest
  julia> Q = Quiver("1--2,1---3,2----3")
  Quiver with adjacency matrix [0 2 3; 0 0 4; 0 0 0]
  ```
  """
  function Quiver(arrows::String)
    pairs = split(arrows, ",")
    pairs = map(p -> split(p, "-"), pairs)
    pairs = map(p -> [parse(Int, p[1]), parse(Int, p[end]), length(p) - 1], pairs)

    n = maximum(maximum(pair[1:2]) for pair in pairs)
    A = zeros(Int, n, n)
    for pair in pairs
      A[pair[1], pair[2]] = pair[3]
    end
    return Quiver(A, "")
  end
end

==(Q1::Quiver, Q2::Quiver) = Q1.adjacency == Q2.adjacency
hash(Q::Quiver) = hash(Q.adjacency)

function show(io::IO, Q::Quiver)
  if Q.name == ""
    print(io, "Quiver with adjacency matrix ")
  else
    print(io, Q.name * ", with adjacency matrix ")
  end
  print(io, Q.adjacency)
end

mutable struct ChowRing
  parent::Any
  ring::Singular.PolyRing{Singular.n_Q}
  chi::AbstractVector{Int}
  point::Union{Singular.spoly{Singular.n_Q},UndefInitializer}
  todd::Union{Singular.spoly{Singular.n_Q},UndefInitializer}
  _R::Singular.PolyRing{Singular.n_Q}
  _inclusion::Singular.SAlgHom{Singular.Rationals}
  ChowRing() = new()
end

function show(io::IO, chow::ChowRing)
  rng = isdefined(chow, :ring) ? chow.ring : UndefInitializer()
  ch = isdefined(chow, :chi) ? chow.chi : UndefInitializer()
  pt = isdefined(chow, :point) ? chow.point : UndefInitializer()
  td = isdefined(chow, :todd) ? chow.todd : UndefInitializer()
  print(
    io,
    "Chow ring on
  $(chow.parent)

  Intersection theory data:

 - Chow ring: $(rng),
 - Linearization: $(ch),
 - Point class: $(pt),
 - Todd class: $(td).
    ",
  )
end

linearization(CH::ChowRing) = CH.chi

"""
# Summary

`abstract type QuiverModuli`

Abstract type for a moduli space or stack of quiver representations.

# Supertype Hierarchy

`QuiverModuliSpace <: QuiverModuli <: Any`\\
`QuiverModuliStack <: QuiverModuli <: Any`

"""
abstract type QuiverModuli end

# TODO consider this:
# https://stackoverflow.com/questions/71738970/in-julia-declare-abstractvectorabstractvector
# this is also necessary to be able to type function outputs correctly.
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
  if condition in ["stable", "semistable"] &&
    length(d) == nvertices(Q) &&
    length(theta) == nvertices(Q)
    d = coerce_vector(d)
    theta = coerce_vector(theta)
    M = QuiverModuliSpace(Q, d, theta, condition, denom, ChowRing())
    setfield!(M.chow, :parent, M)
    return M
  end
  throw(DomainError("Invalid input"))
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
    if condition in ["stable", "semistable"] &&
      length(d) == nvertices(Q) &&
      length(theta) == nvertices(Q)
      d = coerce_vector(d)
      theta = coerce_vector(theta)
      return new(Q, d, theta, condition, denom)
    end
    throw(DomainError("Invalid input"))
  end
end

function show(io::IO, M::QuiverModuliSpace)
  print(
    io,
    "Quiver moduli space defined as follows:
 - quiver: $(M.Q),
 - dimension vector: $(M.d),
 - stability parameter $(M.theta),
 - condition: $(M.condition).
    ",
  )
end
function show(io::IO, M::QuiverModuliStack)
  print(
    io,
    "Quiver moduli stack defined as follows:
 - quiver: $(M.Q),
 - dimension vector: $(M.d),
 - stability parameter $(M.theta),
 - condition: $(M.condition).
    ",
  )
end

"""
# Summary

`struct HNType`

A struct for a Harder-Narasimhan type.

# Fields

 `hn :: Vector{SVector{T,Int}}`\\

"""
struct HNType{T}
  hn::Vector{SVector{T,Int}}
  # should this contain Q, d and slope?
  function HNType(dstar::Vector{<:AbstractVector{Int}})
    T = length(dstar[1])
    return new{T}(coerce_vector.(dstar))
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

"""
# Summary

`mutable struct Bundle`

An abstract bundle on a quiver moduli. It is represented in practice by
its Chern character.

# Fields

`parent :: ChowRing`\\
`rank   :: Int`\\
`chern  :: Singular.spoly{Singular.n_Q}`


"""
mutable struct Bundle
  parent::ChowRing
  rank::Int
  chern_character::Union{Singular.spoly{Singular.n_Q},UndefInitializer}
  chern_class::Union{Dict{Int,Singular.spoly{Singular.n_Q}},UndefInitializer}
  # teleman_weights::Dict{Vector{Any}, Union{Int, Vector{Int}}} # TODO implement
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
  newbundle = Bundle()
  setfield!(newbundle, :parent, parent)
  setfield!(newbundle, :rank, rank)
  cl = Dict{Int,Singular.spoly{Singular.n_Q}}(i => chern_classes[i + 1] for i in 0:n)
  setfield!(newbundle, :chern_class, cl)
  return newbundle
end

function Bundle(parent::ChowRing, rank::Int, chern_class::Singular.spoly{Singular.n_Q})
  n = dimension(parent.parent)
  newbundle = Bundle()
  setfield!(newbundle, :parent, parent)
  setfield!(newbundle, :rank, rank)
  hom = homogeneous_components(parent.parent, chern_class)
  cl = Dict{Int,Singular.spoly{Singular.n_Q}}(i => hom[i + 1] for i in 0:n)
  setfield!(newbundle, :chern_class, cl)
  return newbundle
end

function Bundle(parent::ChowRing, chern_character::Singular.spoly{Singular.n_Q})
  r = constant_coefficient(chern_character)
  denominator(r) != 1 && throw(DomainError("Incorrect Chern character."))
  newbundle = Bundle()
  setfield!(newbundle, :parent, parent)
  setfield!(newbundle, :rank, Int(Singular.numerator(r)))
  setfield!(newbundle, :chern_character, chern_character)
  return newbundle
end
function Bundle(parent::ChowRing, char::Int)
  CH = parent.ring
  newbundle = Bundle()
  setfield!(newbundle, :parent, parent)
  setfield!(newbundle, :rank, char)
  setfield!(newbundle, :chern_character, CH(char))
  return newbundle
end

function Bundle(M::QuiverModuliSpace, char::Singular.spoly{Singular.n_Q})
  newbundle = Bundle()
  setfield!(newbundle, :parent, M.chow)
  r = constant_coefficient(char)
  denominator(r) != 1 && throw(DomainError("Incorrect Chern character."))
  setfield!(newbundle, :rank, Int(Singular.numerator(r)))
  setfield!(newbundle, :chern_character, char)
  return newbundle
end

function Bundle(M::QuiverModuliSpace, rank::Int, x::Singular.spoly{Singular.n_Q})
  newbundle = Bundle()
  setfield!(newbundle, :parent, M.chow)
  setfield!(newbundle, :rank, rank)
  hom = homogeneous_components(M, x)
  cl = Dict{Int,Singular.spoly{Singular.n_Q}}(i => hom[i + 1] for i in 0:dimension(M))
  setfield!(newbundle, :chern_class, cl)
  return newbundle
end

function Bundle(M::QuiverModuliSpace, rank::Int, x::Dict{Int,Singular.spoly{Singular.n_Q}})
  newbundle = Bundle()
  setfield!(newbundle, :parent, M.chow)
  setfield!(newbundle, :rank, rank)
  setfield!(newbundle, :chern_class, x)
  return newbundle
end

function Bundle(M::QuiverModuliSpace, rank::Int, x::Vector{Singular.spoly{Singular.n_Q}})
  newbundle = Bundle()
  setfield!(newbundle, :parent, M.chow)
  setfield!(newbundle, :rank, rank)
  hom = homogeneous_components(M, sum(x))
  cl = Dict{Int,Singular.spoly{Singular.n_Q}}(i => hom[i + 1] for i in 0:dimension(M))
  setfield!(newbundle, :chern_class, x)
  return newbundle
end
