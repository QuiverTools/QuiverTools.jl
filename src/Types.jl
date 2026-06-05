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
  my quiver, with adjacency matrix [0 1; 2 0]
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
      Quiver(arrows)

  Construct a quiver from its arrows encoded in a string.

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
    print(io, Q.adjacency)
  else
    print(io, Q.name)
  end
end

# TODO ChowRing.parent should be of type QuiverModuliSpace,
# but this creates a circular dependency.
"""
# Summary

`struct ChowRing`

A Type used to encode various Chow ring data.

# Fields

`parent :: Any`\\
`ring   :: Singular.PolyRing{Singular.n_Q}`\\
`chi    :: AbstractVector{Int}`\\
`point  :: Union{Singular.spoly{Singular.n_Q},UndefInitializer}`\\
`todd   :: Union{Singular.spoly{Singular.n_Q},UndefInitializer}`\\
`_R     :: Singular.PolyRing{Singular.n_Q}`\\
`_inclusion :: Singular.SAlgHom{Singular.Rationals}`
"""
mutable struct ChowRing
  parent::Any
  ring::Singular.PolyRing{Singular.n_Q}
  chi::AbstractVector{Int}
  point::Union{Singular.spoly{Singular.n_Q},UndefInitializer}
  todd::Union{Singular.spoly{Singular.n_Q},UndefInitializer}
  _R::Singular.PolyRing{Singular.n_Q}
  _inclusion::Singular.SAlgHom{Singular.Rationals}

  """   ChowRing()"""
  ChowRing() = new()
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
    length(d) == n_vertices(Q) &&
    length(theta) == n_vertices(Q)
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
      length(d) == n_vertices(Q) &&
      length(theta) == n_vertices(Q)
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

# TODO this needs to be explained better
"""
# Summary

`struct HNType`

A struct for a Harder-Narasimhan type.

# Fields

 `hn :: Vector{AbstractVector{Int}}`\\

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

A struct to encode Luna types.

# Fields

 `data :: Dict{AbstractVector{Int},AbstractVector{Int}}`\\

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
