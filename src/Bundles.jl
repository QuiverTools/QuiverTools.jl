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
  chern_character::Singular.spoly{Singular.n_Q}
  chern_class::Singular.spoly{Singular.n_Q}
  # teleman_weights::Dict{Vector{Any}, Union{Int, Vector{Int}}} # TODO implement
end

function show(io::IO, F::Bundle)
  print(
    io,
    "Bundle of rank $(rank(F)), with Chern character
$(F.chern_character)",
  )
end

function Bundle(parent::ChowRing, chern_character::Singular.spoly{Singular.n_Q})
  r = Int(numerator(QuiverTools.constant_coefficient(chern_character)))
  return Bundle(parent, r, chern_character)
end
Bundle(parent::ChowRing, char::Int) = Bundle(parent, char, parent.ring(char))
Bundle(M::QuiverModuliSpace, char) = Bundle(M.chow, char)

rank(F::Bundle) = F.rank
chern_character(F::Bundle) = F.chern_character

function chern_classes(F::Bundle)
  !isdefined(F, :chern_class) && setfield!(F, :chern_class, _chern_classes_from_character(F))
  return F.chern_class
end
function chern_class(F::Bundle)
  !isdefined(F, :chern_class) && setfield!(F, :chern_class, _chern_classes_from_character(F))
  return sum(chern_classes(F))
end

function chern_class(F::Bundle, k)
  !isdefined(F, :chern_class) && setfield!(F, :chern_class, _chern_classes_from_character(F))
  return chern_classes(F)[k+1]
end
chow_ring(F::Bundle) = F.parent.ring
variety(F::Bundle) = F.parent.parent
structure_sheaf(M::QuiverModuliSpace) = Bundle(M, 1)

##############################
# Operations on Bundle objects
##############################

"""
    dual(F::Bundle)

Return the dual bundle of `F`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M))
Bundle of rank 1, with Chern character
x21 + 1

julia> dual(td)
Bundle of rank 1, with Chern character
-x21 + 1
```

On our favourite 6-fold:

```jldoctest
julia> Q, d = kronecker_quiver(3), [2, 3];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=[-1, 1]);

julia> td = Bundle(M, todd_class(M))
Bundle of rank 1, with Chern character
-17//8*x12*x21 + x21^2 + 823//360*x12*x22 - 823//1080*x22^2 + 553//1080*x21*x23 - 77//60*x22*x23 + x23^2 + 5//12*x12 - 3//2*x21 + 9//8*x23 + 1


julia> dual(td)
Bundle of rank 1, with Chern character
17//8*x12*x21 + x21^2 + 823//360*x12*x22 - 823//1080*x22^2 + 553//1080*x21*x23 + 77//60*x22*x23 + x23^2 + 5//12*x12 + 3//2*x21 - 9//8*x23 + 1
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

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M))
Bundle of rank 1, with Chern character
x21 + 1

julia> F = 3*td
Bundle of rank 3, with Chern character
3*x21 + 3

julia> map(i -> exterior_power(F, i), 0:4)
5-element Vector{Bundle}:
 Bundle of rank 1, with Chern character
1
 Bundle of rank 3, with Chern character
3*x21 + 3
 Bundle of rank 3, with Chern character
6*x21 + 3
 Bundle of rank 1, with Chern character
3*x21 + 1
 Bundle of rank 0, with Chern character
0
```
"""
function exterior_power(F::Bundle, k::Int)
  return Bundle(F.parent, _chern_characters_wedge(F, k)[end])
end
det(F::Bundle) = exterior_power(F, rank(F))

"""
    symmetric_power(F::Bundle, k::Int)

Return the `k`-th symmetric power of `F`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M)); F = 3*td
Bundle of rank 3, with Chern character
3*x21 + 3

julia> map(i -> symmetric_power(F, i), 0:4)
5-element Vector{Bundle}:
 Bundle of rank 1, with Chern character
1
 Bundle of rank 3, with Chern character
3*x21 + 3
 Bundle of rank 6, with Chern character
12*x21 + 6
 Bundle of rank 10, with Chern character
30*x21 + 10
 Bundle of rank 15, with Chern character
60*x21 + 15
```
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
      CH(1//j) * truncate(M,
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
  r = rank(F)

  wedges = _chern_characters_wedge(F, r)
  # init as CH(0) for type stability
  syms = [CH(0) for _ in 1:(k + 1)]
  syms[1], syms[2] = CH(1), x

  for j in 2:k
    syms[j + 1] = truncate(M,
      sum(
        (-CH(1))^(i + 1) * wedges[i + 1] * syms[j - i + 1] for i in 1:min(j, r);
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


function _chern_classes_from_character(F::Bundle)
  CH = chow_ring(F)
  M = variety(F)
  n = dimension(M)
  comps = homogeneous_components(M, chern_character(F))
  p = [(CH(-1))^i * factorial(CH(i)) * comps[i+1] for i in 0:n]
  e = [CH(0) for _ in n+1]
  e[1] = CH(1)
  for i in 1:n
    e[i+1] = CH(-1//i) * sum(p[j+1] * e[i-j+1] for j in 1:i)
    # e[i+1] = div(e[i+1], CH(1)) # simplify
  end
  return e
end
