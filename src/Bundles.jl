##################################
# Tensor calculus on quiver moduli
##################################

_has_chern_data(F::Bundle) = isdefined(F, :chern_character) || isdefined(F, :chern_class)
"""   chern_character(F::Bundle)"""
function chern_character(F::Bundle)
  !_has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_character) &&
    setfield!(F, :chern_character, _chern_character_from_classes(F))
  return F.chern_character
end
"""   chern_classes(F::Bundle)"""
function chern_classes(F::Bundle)
  !_has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, _chern_classes_from_character(F))
  return F.chern_class
end
"""   chern_class(F::Bundle)"""
function chern_class(F::Bundle)
  !_has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, _chern_classes_from_character(F))
  return sum(values(chern_classes(F)))
end
"""   chern_class(F::Bundle, k)"""
function chern_class(F::Bundle, k)
  !_has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, _chern_classes_from_character(F))
  return chern_classes(F)[k]
end
"""    teleman_weights(F::Bundle)"""
function teleman_weights(F::Bundle)
  !isdefined(F, :teleman_weights) && throw(ArgumentError("Bundle has no weights."))
  return F.teleman_weights
end
"""    rank(F::Bundle)"""
rank(F::Bundle) = F.rank

"""   chow_ring(F::Bundle)"""
chow_ring(F::Bundle) = F.parent.ring

"""    variety(F::Bundle)"""
variety(F::Bundle) = F.parent.parent

"""    structure_sheaf(M::QuiverModuliSpace)"""
function structure_sheaf(M::QuiverModuliSpace)
  new = Bundle(M, 1)
  HN = all_hn_types(M; unstable=true, ordered=false)
  set_teleman_weights!(new, Dict(hn_type => [0] for hn_type in HN))
  return new
end

"""    zero_sheaf(M::QuiverModuliSpace)"""
function zero_sheaf(M::QuiverModuliSpace)
  new = Bundle(M, 0)
  HN = all_hn_types(M; unstable=true, ordered=false)
  set_teleman_weights!(new, Dict(hn_type => [] for hn_type in HN))
  return new
end

##############################
# Operations on Bundle objects
##############################

"""
    dual(F::Bundle)

Compute the dual bundle of `F`.

# Input

- `F::Bundle`: a bundle.

# Output

- the dual bundle of `F`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M))
Bundle of rank 1

julia> chern_character(dual(td))
-x21 + 1
```

On our favourite 6-fold:

```jldoctest
julia> Q, d = kronecker_quiver(3), [2, 3];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=[-1, 1]);

julia> td = Bundle(M, todd_class(M));

julia> chern_character(td)
-17//8*x12*x21 + x21^2 + 823//360*x12*x22 - 823//1080*x22^2 + 553//1080*x21*x23 - 77//60*x22*x23 + x23^2 + 5//12*x12 - 3//2*x21 + 9//8*x23 + 1


julia> tdd = dual(td)
Bundle of rank 1

julia> chern_character(tdd)
17//8*x12*x21 + x21^2 + 823//360*x12*x22 - 823//1080*x22^2 + 553//1080*x21*x23 + 77//60*x22*x23 + x23^2 + 5//12*x12 + 3//2*x21 - 9//8*x23 + 1
```
"""
function dual(F::Bundle)
  if _has_chern_data(F)
    new = Bundle(F.parent, adams(F, -1))
  else
    new = Bundle()
    setfield!(new, :parent, F.parent)
    setfield!(new, :rank, F.rank)
  end
  if isdefined(F, :teleman_weights)
    weights_dual = Dict(
      hn_type => -teleman_weights(F)[hn_type] for hn_type in keys(teleman_weights(F))
    )
    set_teleman_weights!(new, weights_dual)
  end
  return new
end

function *(n::Int, F::Bundle)
  n == 0 && return zero_sheaf(variety(F))
  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, F.rank * n)
  _has_chern_data(F) && setfield!(new, :chern_character, n * chern_character(F))
  if isdefined(F, :teleman_weights)
    n_weights = Dict(
      hn_type => reduce(vcat, teleman_weights(F)[hn_type] for _ in 1:n)
      for hn_type in keys(teleman_weights(F))
    )
    set_teleman_weights!(new, n_weights)
  end
  return new
end

*(F::Bundle, n::Int) = n * F

# power = tensor product
function ^(F::Bundle, n::Int)
  n == 0 && return structure_sheaf(variety(F))
  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, (F.rank)^n)
  _has_chern_data(F) && setfield!(new, :chern_character, chern_character(F)^n)
  if isdefined(F, :teleman_weights)
    pow_weights = Dict(
      hn_type =>
        map(sum, Iterators.product(map(k -> teleman_weights(F)[hn_type], 1:n)...))[:]
      for hn_type in keys(teleman_weights(F))
    )
    set_teleman_weights!(new, pow_weights)
  end
  return new
end

function +(F::Bundle, G::Bundle)
  F.parent != G.parent && throw(DomainError("Different Chow rings."))
  homog_chow = _has_chern_data(F) == _has_chern_data(G)
  !homog_chow && throw(ArgumentError("Dishomogeneous Chern data."))
  homog_weights = isdefined(F, :teleman_weights) == isdefined(G, :teleman_weights)
  !homog_weights && throw(ArgumentError("Dishomogeneous Teleman weights."))

  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, F.rank + G.rank)

  _has_chern_data(F) &&
    setfield!(new, :chern_character, chern_character(F) + chern_character(G))
  if isdefined(F, :teleman_weights)
    new_weights = Dict(
      hn_type => vcat(teleman_weights(F)[hn_type], teleman_weights(G)[hn_type])
      for hn_type in keys(teleman_weights(F))
    )
    set_teleman_weights!(new, new_weights)
  end
  return new
end

function -(F::Bundle, G::Bundle)
  F.parent != G.parent && throw(DomainError("Different Chow rings."))
  (isdefined(F, :teleman_weights) || isdefined(G, :teleman_weights)) && throw(
    ArgumentError("Cannot subtract bundles with weights.")
  )
  Bundle(F.parent, chern_character(F) - chern_character(G))
end

function *(F::Bundle, G::Bundle)
  F.parent != G.parent && throw(DomainError("Different Chow rings."))
  homog_chow = _has_chern_data(F) == _has_chern_data(G)
  !homog_chow && throw(ArgumentError("Dishomogeneous Chern data."))
  homog_weights = isdefined(F, :teleman_weights) == isdefined(G, :teleman_weights)
  !homog_weights && throw(ArgumentError("Dishomogeneous Teleman weights."))

  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, F.rank * G.rank)

  _has_chern_data(F) &&
    setfield!(new, :chern_character, chern_character(F) * chern_character(G))
  if isdefined(F, :teleman_weights)
    new_weights = Dict(
      hn_type =>
        [x + y for x in teleman_weights(F)[hn_type] for y in teleman_weights(G)[hn_type]]
      for hn_type in keys(teleman_weights(F))
    )
    set_teleman_weights!(new, new_weights)
  end
  return new
end

"""
    exterior_power(F::Bundle, k::Int)

Compute the `k`-th exterior power of `F`.

# Input

- `F::Bundle`: a bundle.
- `k::Int`: the degree of the exterior power.

# Output

- the `k`-th exterior power of `F`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M))
Bundle of rank 1

julia> F = 3*td
Bundle of rank 3

julia> chern_character(F)
3*x21 + 3

julia> W = map(i -> exterior_power(F, i), 0:4);

julia> map(w -> (rank(w), chern_character(w)), W)
5-element Vector{Tuple{Int64, Singular.spoly{Singular.n_Q}}}:
 (1, 1)
 (3, 3*x21 + 3)
 (3, 6*x21 + 3)
 (1, 3*x21 + 1)
 (0, 0)
```
"""
function exterior_power(F::Bundle, k::Int)
  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, binomial(rank(F), k))
  CH = chow_ring(F)

  _has_chern_data(F) &&
    setfield!(new, :chern_character, simplify!(CH(_chern_characters_wedge(F, k)[end])))
  if isdefined(F, :teleman_weights)
    new_weights = Dict(
      hn_type =>
        [sum(c) for c in combinations(teleman_weights(F)[hn_type], k)]
      for hn_type in keys(teleman_weights(F))
    )
    set_teleman_weights!(new, new_weights)
  end
  return new
end

"""
    det(F::Bundle)

Compute the determinant of `F`. This is the top exterior power of `F`.
"""
det(F::Bundle) = exterior_power(F, rank(F))

"""
    symmetric_power(F::Bundle, k::Int)

Compute the `k`-th symmetric power of `F`.

# Input

- `F::Bundle`: a bundle.
- `k::Int`: the degree of the symmetric power.

# Output

- the `k`-th symmetric power of `F`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> td = Bundle(M, todd_class(M));

julia> F = 3*td
Bundle of rank 3

julia> chern_character(F)
3*x21 + 3

julia> W = map(i -> symmetric_power(F, i), 0:4);

julia> map(w -> (rank(w), chern_character(w)), W)
5-element Vector{Tuple{Int64, Singular.spoly{Singular.n_Q}}}:
 (1, 1)
 (3, 3*x21 + 3)
 (6, 12*x21 + 6)
 (10, 30*x21 + 10)
 (15, 60*x21 + 15)
```
"""
function symmetric_power(F::Bundle, k::Int)
  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, binomial(rank(F) + k - 1, rank(F) - 1))
  CH = chow_ring(F)

  _has_chern_data(F) &&
    setfield!(new, :chern_character, simplify!(CH(_chern_characters_symmetric(F, k)[end])))
  if isdefined(F, :teleman_weights)
    new_weights = Dict(
      hn_type =>
        [sum(c) for c in with_replacement_combinations(teleman_weights(F)[hn_type], k)]
      for hn_type in keys(teleman_weights(F))
    )
    set_teleman_weights!(new, new_weights)
  end
  return new
end

function homogeneous_components(M::QuiverModuliSpace, x; unsafe::Bool=false)
  if unsafe
    n = 1 - euler_form(M.Q, M.d, M.d)
  else
    n = dimension(M)
  end

  CH = chow_ring(M; unsafe=unsafe)
  return [
    sum(
      t for t in Singular.terms(x) if __chow_ring_monomial_grading(M, t) == i; init=CH(0)
    )
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
    simplify!(wedges[j + 1])
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
    simplify!(syms[j + 1])
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
  p = [(CH(-1))^i * CH(factorial(i)) * comps[i + 1] for i in 0:n]
  e = [CH(0) for _ in 1:(n + 1)]
  e[1] = CH(1)
  for i in 1:n
    e[i + 1] = CH(-1//i) * sum(p[j + 1] * e[i - j + 1] for j in 1:i)
    simplify!(e[i + 1])
  end
  return Dict(i => e[i + 1] for i in 0:n)
end

function _chern_character_from_classes(F::Bundle)
  CH = F.parent.ring
  M = variety(F)
  n = dimension(M)
  n == 0 && return CH(0)
  e = chern_classes(F)
  p = vcat([-e[1]], [CH(0) for _ in 1:(n - 1)])
  for i in 1:(n - 1)
    p[i + 1] = -CH(i + 1) * e[i + 1] - sum(e[j] * p[i - j + 1] for j in 1:i)
  end
  return simplify(sum(CH((-1)^i//factorial(i)) * p[i] for i in 1:n) + rank(F))
end

simplify(f::Singular.spoly{Singular.n_Q}) = div(f, f.parent(1))

function simplify!(f::Singular.spoly{Singular.n_Q})
  f = div(f, f.parent(1))
  return f
end

"""
    line_bundle(M::QuiverModuliSpace, eta::AbstractVector{Int}; teleman::Bool=true)

Construct the descent of `L(eta)` on the quiver moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a quiver moduli space.
- `eta::AbstractVector{Int}`: a dimension vector.

# Output

- the line bundle on `M` with Chern class `chern_class_line_bundle(M, eta)`.

# Examples

The ample generator of the Picard group of our favourite 6-fold:

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> chow_ring(M; chi=[-1, 1]);

julia> H = line_bundle(M, [3, -2]);

julia> chern_class(H)
-x21

julia> teleman_weights(H)
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [-10]
  [[2, 1], [0, 2]]         => [-20]
  [[1, 0], [1, 2], [0, 1]] => [-40]
  [[1, 0], [1, 3]]         => [-15]
  [[1, 0], [1, 1], [0, 2]] => [-35]
  [[1, 1], [1, 2]]         => [-5]
  [[2, 0], [0, 3]]         => [-30]
```
"""
function line_bundle(M::QuiverModuliSpace, eta::AbstractVector{Int}; teleman::Bool=true)
  eta' * M.d != 0 && throw(ArgumentError("$(eta) is not a linearization."))
  new = Bundle(M, 1, chern_class_line_bundle(M, eta))
  teleman && set_teleman_weights!(new, weights_line_bundle(M, eta))
  return new
end

"""
    canonical_bundle(M::QuiverModuliSpace; teleman::Bool=true, verbose::Bool=false, unsafe::Bool=false)

Compute the canonical bundle on the quiver moduli space `M`.

If the moduli problem is amply stable and does not admit properly semistable representations,
the canonical bundle is described in
[[Proposition 4.2, MR4352662](https://mathscinet.ams.org/mathscinet-getitem?mr=4352662)].

This function computes both the Chern character and the Teleman weights
of the canonical bundle.

# Input

- `M::QuiverModuliSpace`: a quiver moduli space.
- `teleman::Bool`: Optional keyword argument to compute the Teleman weights of the canonical bundle. Default is `true`.
- `verbose::Bool`: Optional keyword argument to print warnings. Default is `false`.
- `unsafe::Bool`: Optional keyword argument to skip ample stability and properly semistable checks. Default is `false`.

# Output

- the canonical bundle on `M`.

# Example

On various projective spaces:

```jldoctest
julia> d = [1, 1]; Pn = map(i -> QuiverModuliSpace(kronecker_quiver(i + 1), d), 1:5);

julia> omega = map(canonical_bundle, Pn);

julia> map(chern_class, omega)
5-element Vector{Singular.spoly{Singular.n_Q}}:
 2*x11
 3*x11
 4*x11
 5*x11
 6*x11

julia> map(teleman_weights, omega)
5-element Vector{Dict{HNType, Vector{Int64}}}:
 Dict([[1, 0], [0, 1]] => [4])
 Dict([[1, 0], [0, 1]] => [6])
 Dict([[1, 0], [0, 1]] => [8])
 Dict([[1, 0], [0, 1]] => [10])
 Dict([[1, 0], [0, 1]] => [12])
```

On our favourite 6-fold:

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> omega = canonical_bundle(M);

julia> chern_class(omega)
3*x21

julia> QuiverTools.teleman_weights(omega)
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [30]
  [[2, 1], [0, 2]]         => [60]
  [[1, 0], [1, 2], [0, 1]] => [120]
  [[1, 0], [1, 3]]         => [45]
  [[1, 0], [1, 1], [0, 2]] => [105]
  [[1, 1], [1, 2]]         => [15]
  [[2, 0], [0, 3]]         => [90]
```
"""
function canonical_bundle(M::QuiverModuliSpace; verbose::Bool=false, unsafe::Bool=false)
  if !unsafe
    has_properly_semistables(M.Q, M.d, M.theta, M.denom) &&
      throw(
        ArgumentError(
          "The quiver moduli problem has properly semistables, no description of the canonical bundle is available."
        ),
      )
    !is_amply_stable(M) &&
      throw(
        ArgumentError(
          "The quiver moduli problem is not amply stable, no description of the canonical bundle is available."
        ),
      )
  else
    verbose && @warn "Unsafe computation."
  end
  return line_bundle(M, -canonical_stability(M.Q, M.d))
end

"""
    universal_bundle(M:::QuiverModuliSpace, i::Int)

Compute the `i`-th universal bundle of `M`.

This function computes both the Chern classes and the Teleman weights
of the universal bundle.

To use this method, the moduli space `M` must have a linearization.
If not, this method will initialize the Chow ring of `M`
by calling `chow_ring(M)` and it will use the default linearization.

# Input

- `M::QuiverModuliSpace`: a quiver moduli space.
- `i::Int`: the universal bundle on the `i`-th vertex of the quiver.
- `teleman::Bool`: Optional keyword argument to compute the Teleman weights of the universal bundle. Default is `true`.
- `unsafe::Bool`: Optional keyword argument to skip ample stability checks. Default is `false`.

# Output

- the `i`-th universal bundle on `M`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=[1, 0]);

julia> u1, u2 = universal_bundle(M, 1), universal_bundle(M, 2);

julia> map(chern_class, [u1, u2])
2-element Vector{Singular.spoly{Singular.n_Q}}:
 x11 + 1
 x21 + 1

julia> map(teleman_weights, [u1, u2])
2-element Vector{Dict{HNType, Vector{Int64}}}:
 Dict([[1, 0], [0, 1]] => [0])
 Dict([[1, 0], [0, 1]] => [-2])
```

On our favourite 6-fold:

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> chow_ring(M; chi=[-1, 1]);

julia> u1, u2 = universal_bundle(M, 1), universal_bundle(M, 2);

julia> map(chern_character, [u1, u2])
2-element Vector{Singular.spoly{Singular.n_Q}}:
 1//6*x12*x21 + 1//2*x21^2 + 1//12*x12*x22 - 1//36*x22^2 - 7//72*x21*x23 - 1//120*x22*x23 - 1//720*x23^2 - x12 + x21 - 1//2*x23 + 2
 2//3*x12*x21 + 1//2*x21^2 - 1//6*x12*x22 - 1//2*x21*x22 + 1//12*x22^2 + 1//24*x21*x23 + 1//180*x22*x23 + x21 - x22 + 3

julia> w = map(teleman_weights, [u1, u2]);

julia> w[1]
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [5, 5]
  [[2, 1], [0, 2]]         => [10, 10]
  [[1, 0], [1, 2], [0, 1]] => [25, 15]
  [[1, 0], [1, 3]]         => [10, 5]
  [[1, 0], [1, 1], [0, 2]] => [20, 15]
  [[1, 1], [1, 2]]         => [5, 0]
  [[2, 0], [0, 3]]         => [15, 15]

julia> w[2]
Dict{HNType, Vector{Int64}} with 7 entries:
  [[2, 2], [0, 1]]         => [5, 5, 0]
  [[2, 1], [0, 2]]         => [10, 5, 5]
  [[1, 0], [1, 2], [0, 1]] => [15, 15, 10]
  [[1, 0], [1, 3]]         => [5, 5, 5]
  [[1, 0], [1, 1], [0, 2]] => [15, 10, 10]
  [[1, 1], [1, 2]]         => [5, 0, 0]
  [[2, 0], [0, 3]]         => [10, 10, 10]
```
"""
function universal_bundle(
  M::QuiverModuliSpace, i::Int; unsafe::Bool=false, teleman::Bool=true
)
  gcd(M.d) > 1 && throw(
    ArgumentError("gcd($(M.d))  = $(gcd(M.d)) > 1, the universal bundles do not exist.")
  )
  cl = total_chern_class_universal(M, i; unsafe=unsafe)
  new = Bundle(M, M.d[i], cl)

  teleman && set_teleman_weights!(new, weights_universal_bundle(M, i))
  return new
end

"""
    degree(F::Bundle; unsafe::Bool=false)

Compute the degree of the bundle `F`.
If `rank(F)` is larger than ``1``, returns the degree of the determinant of `F`.

# Input

- `F::Bundle`: a bundle.
- `unsafe::Bool=false`: whether to skip ample stability checks. Default is `false`.

# Output

- the degree of `F`.

# Example

The degrees of canonical bundles on the first projective spaces:

```jldoctest
julia> d = [1, 1]; Pn = map(i -> QuiverModuliSpace(kronecker_quiver(i + 1), d), 1:5);

julia> omega = map(canonical_bundle, Pn);

julia> omega = map(dual, omega);

julia> map(degree, omega)
5-element Vector{Singular.spoly{Singular.n_Q}}:
 2
 9
 64
 625
 7776
```

An example from [[arXiv:2411.15125](https://arxiv.org/abs/2411.15125)]:

```jldoctest
julia> Q = Quiver("1-2,1-3,2---3"); d = [1, 1, 1]; a = [1, 1, -1];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> F = dual(canonical_bundle(M))
Bundle of rank 1

julia> chern_class(F)
2*x31 + 1

julia> degree(F)
56
```
"""
function degree(F::Bundle; unsafe::Bool=false)
  M = variety(F)
  n = dimension(M)

  out = chern_class(det(F))
  out = Singular.jet(out, n)

  if n >= 2
    class_det = deepcopy(out)
    for _ in 1:(n - 1)
      # the multiplication here takes most of runtime
      Oscar.mul!(out, out, class_det)
      out = Singular.jet(out, n)
    end
  end

  pt = point_class(M; unsafe=unsafe)
  return div(homogeneous_components(M, out)[n + 1], pt)
end
