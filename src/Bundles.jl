###############################################################################
# Tensor calculus on quiver moduli
###############################################################################

import Base: *, +, -, ^

export Bundle, chern_character, dual, exterior_power, symmetric_power, det

"""
# Summary

`struct Bundle`

An abstract bundle on a quiver moduli. It is represented in practice by
its Chern character.

# Fields

`parent :: ChowRing`\\
`rank   :: Int`\\
`chern  :: Singular.spoly{Singular.n_Q}`
"""
struct Bundle
  parent::ChowRing
  rank::Int
  chern::Singular.spoly{Singular.n_Q}
  # teleman_weights::Dict{Vector{Any}, Union{Int, Vector{Int}}} # TODO implement
end

function show(io::IO, F::Bundle)
  print(
    io,
    "Bundle of rank $(F.rank), with Chern character
$(F.chern)",
  )
end

function Bundle(parent::ChowRing, chern::Singular.spoly{Singular.n_Q})
  r = Int(numerator(QuiverTools.constant_coefficient(chern)))
  return Bundle(parent, r, chern)
end
Bundle(parent::ChowRing, chern::Int) = Bundle(parent, chern, parent.ring(chern))
Bundle(M::QuiverModuliSpace, chern) = Bundle(M.chow, chern)

chern_character(F::Bundle) = F.chern
chow_ring(F::Bundle) = F.parent.ring
variety(F::Bundle) = F.parent.parent

##############################
# Operations on Bundle objects
##############################

"""
    dual(F::Bundle)

Return the dual bundle of `F`.

```
"""
function dual(F::Bundle)
  return Bundle(F.parent, adams(F, -1))
end

*(n::Int, F::Bundle) = Bundle(F.parent, n * chern_character(F))
*(F::Bundle, n::Int) = n * F
^(F::Bundle, n::Int) = Bundle(F.parent, chern_character(F)^n)

# direct sum, quotient and tensor product
+(F::Bundle, G::Bundle) =
  if F.parent == G.parent
    Bundle(F.parent, chern_character(F) + chern_character(G))
  else
    throw(DomainError("Different Chow rings."))
  end
-(F::Bundle, G::Bundle) =
  if F.parent == G.parent
    Bundle(F.parent, chern_character(F) - chern_character(G))
  else
    throw(DomainError("Different Chow rings."))
  end
*(F::Bundle, G::Bundle) =
  if F.parent == G.parent
    Bundle(F.parent, chern_character(F) * chern_character(G))
  else
    throw(DomainError("Different Chow rings."))
  end

"""
    exterior_power(F::Bundle, k::Int)

Return the `k`-th exterior power of `F`.
"""
function exterior_power(F::Bundle, k::Int)
  return Bundle(F.parent, _chern_characters_wedge(F, k)[end])
end
det(F::Bundle) = exterior_power(F, F.rank)

"""
    symmetric_power(F::Bundle, k::Int)

Return the `k`-th symmetric power of `F`.
"""
function symmetric_power(F::Bundle, k::Int)
  return Bundle(F.parent, _chern_characters_symmetric(F, k)[end])
end

function homogeneous_components(M::QuiverModuliSpace, x)
  n = dimension(M)
  return [
    sum(t for t in Singular.terms(x) if __chow_ring_monomial_grading(M, t) == i; init=0)
    for
    i in 0:n
  ]
end

function truncate(M::QuiverModuliSpace, x, n)
  comps = homogeneous_components(M, x)
  # TODO this may be slow, try building new polynomial
  # also can I use the incorrect but faster Singular.degree?
  return sum(comps[i + 1] for i in 0:n)
end

"""
    _chern_characters_wedge(F::Bundle, k)

Compute the exterior powers of `F` up to degree `k`.
For internal use only.
"""
function _chern_characters_wedge(F::Bundle, k)
  k == 0 && return [1]
  x = chern_character(F)
  M = variety(F)
  CH = chow_ring(F)
  n = dimension(M)

  # init as CH(0) for type stability
  wedges = [CH(0) for _ in 1:(k + 1)]
  wedges[1], wedges[2] = CH(1), x

  for j in 2:k
    wedges[j + 1] =
      1//j * truncate(M,
        sum(
          (-CH(1))^(j - i + 1) * wedges[i + 1] * adams(F, j - i) for i in 0:(j - 1);
          init=CH(0),
        ),
        n)
  end
  return wedges
end

"""
    _chern_characters_symmetric(F::Bundle, k)

Compute the symmetric powers of `F` up to degree `k`.
For internal use only.
"""
function _chern_characters_symmetric(F::Bundle, k)
  k == 0 && return [1]
  x = chern_character(F)
  M = variety(F)
  n = dimension(M)
  CH = chow_ring(F)
  r = F.rank

  wedges = _chern_characters_wedge(F, r)
  # init as CH(0) for type stability
  syms = [CH(0) for _ in 1:(k + 1)]
  syms[1], syms[2] = CH(1), x

  for j in 2:k
    syms[j + 1] = truncate(M,
      sum(
        (-CH(1))^(i + 1) * wedges[i + 1] * syms[j - i + 1] for i in 1:minimum(j, r);
        init=CH(0),
      ),
      n)
  end
  return syms
end

"""
    adams(F::Bundle, k)

Compute the Adams operation ``\\Phi^k`` on the Chern character of `F`.
For internal use only.
"""
function adams(F::Bundle, k)
  M = variety(F)
  n = dimension(M)
  x = chern_character(F)

  return [k^i for i in 0:n]' * homogeneous_components(M, x)
end
