##################################
# Tensor calculus on quiver moduli
##################################

__has_chern_data(F::Bundle) = isdefined(F, :chern_character) || isdefined(F, :chern_class)
"""   chern_character(F::Bundle)"""
function chern_character(F::Bundle)
  !__has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_character) &&
    setfield!(F, :chern_character, __chern_character_from_classes(F))
  return F.chern_character
end
"""   chern_classes(F::Bundle)"""
function chern_classes(F::Bundle)
  !__has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, __chern_classes_from_character(F))
  return F.chern_class
end
"""   chern_class(F::Bundle)"""
function chern_class(F::Bundle)
  !__has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, __chern_classes_from_character(F))
  return sum(values(chern_classes(F)))
end
"""   chern_class(F::Bundle, k)"""
function chern_class(F::Bundle, k)
  !__has_chern_data(F) &&
    throw(
      ArgumentError(
        "Bundle has no Chow ring data."
      ),
    )
  !isdefined(F, :chern_class) &&
    setfield!(F, :chern_class, __chern_classes_from_character(F))
  return chern_classes(F)[k]
end

"""
    segre_class(F::Bundle)

Compute the total Segre class of the bundle `F`, i.e. the inverse of its total Chern
class in the Chow ring, from [`chern_class(F)`](@ref chern_class).

# Examples

The tangent bundle of the projective plane has total Chern class ``1 + 3H + 3H^2``,
hence Segre class ``1 - 3H + 6H^2``. Here the hyperplane class is ``H = -x_{1,1}``:

```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [1, 1]);

julia> T = tangent_bundle(M);

julia> chern_class(T)
3*x1_1^2 - 3*x1_1 + 1

julia> segre_class(T)
6*x1_1^2 + 3*x1_1 + 1
```
"""
segre_class(F::Bundle) = segre_class(variety(F), chern_class(F))

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

"""
    structure_sheaf(M::QuiverModuliSpace)

Construct the structure sheaf of `M`, i.e. the trivial line bundle, with Teleman weight
``0`` on every HN stratum.

# Example

```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> O = structure_sheaf(M);

julia> rank(O), chern_character(O)
(1, 1)

julia> teleman_weights(O) ==
        Dict{HNType, Vector{Int64}}(
             [[2, 2], [0, 1]]         => [0],
             [[2, 1], [0, 2]]         => [0],
             [[1, 0], [1, 2], [0, 1]] => [0],
             [[1, 0], [1, 3]]         => [0],
             [[1, 0], [1, 1], [0, 2]] => [0],
             [[1, 1], [1, 2]]         => [0],
             [[2, 0], [0, 3]]         => [0]
       )
true
```
"""
function structure_sheaf(M::QuiverModuliSpace)
  chow_ring(M)
  new = Bundle(M, 1)
  HN = all_hn_types(M; unstable=true, ordered=false)
  set_teleman_weights!(new, Dict(hn_type => [0] for hn_type in HN))
  return new
end

"""
    zero_sheaf(M::QuiverModuliSpace)

Construct the zero sheaf on `M`, of rank ``0`` and without Teleman weights.

# Example

```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> Z = QuiverTools.zero_sheaf(M);

julia> rank(Z), chern_character(Z)
(0, 0)

julia> teleman_weights(Z) ==
        Dict{HNType, Vector{Int64}}(
             [[2, 2], [0, 1]]         => [],
             [[2, 1], [0, 2]]         => [],
             [[1, 0], [1, 2], [0, 1]] => [],
             [[1, 0], [1, 3]]         => [],
             [[1, 0], [1, 1], [0, 2]] => [],
             [[1, 1], [1, 2]]         => [],
             [[2, 0], [0, 3]]         => []
       )
true
```
"""
function zero_sheaf(M::QuiverModuliSpace)
  chow_ring(M)
  new = Bundle(M, 0)
  HN = all_hn_types(M; unstable=true, ordered=false)
  set_teleman_weights!(new, Dict(hn_type => Int[] for hn_type in HN))
  return new
end

"""
    tangent_bundle(M::QuiverModuliSpace; unsafe::Bool=false)

Compute the tangent bundle of the quiver moduli space `M`.

Uses the identity
``[\\mathrm{T}_M] = \\sum_{a: i \\to j} [U_i^\\vee \\otimes U_j] - \\sum_i [U_i^\\vee \\otimes U_i] + [\\mathcal{O}_M]``
in the Grothendieck group of `M`.

# Example

The tangent bundle of `P^3`, realized as the moduli space of the
generalized Kronecker quiver with 4 arrows and dimension vector `(1, 1)`:

```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(4), [1, 1]); chow_ring(M);

julia> T = tangent_bundle(M)
Bundle of rank 3

julia> chern_class(T, 1)
-4*x1_1
```
"""
tangent_bundle(M::QuiverModuliSpace; unsafe::Bool=false) =
  Bundle(M, chern_character_tangent_bundle(M; unsafe=unsafe))

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
-x2_1 + 1
```

On our favourite 6-fold:

```jldoctest
julia> Q, d = kronecker_quiver(3), [2, 3];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=[-1, 1]);

julia> td = Bundle(M, todd_class(M));

julia> chern_character(td)
-17//8*x1_2*x2_1 + x2_1^2 + 823//360*x1_2*x2_2 - 823//1080*x2_2^2 + 553//1080*x2_1*x2_3 - 77//60*x2_2*x2_3 + x2_3^2 + 5//12*x1_2 - 3//2*x2_1 + 9//8*x2_3 + 1


julia> tdd = dual(td)
Bundle of rank 1

julia> chern_character(tdd)
17//8*x1_2*x2_1 + x2_1^2 + 823//360*x1_2*x2_2 - 823//1080*x2_2^2 + 553//1080*x2_1*x2_3 + 77//60*x2_2*x2_3 + x2_3^2 + 5//12*x1_2 + 3//2*x2_1 - 9//8*x2_3 + 1
```
"""
function dual(F::Bundle)
  if __has_chern_data(F)
    new = Bundle(F.parent, dual(variety(F), chern_character(F)))
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
  __has_chern_data(F) && setfield!(new, :chern_character, n * chern_character(F))
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
  if __has_chern_data(F)
    # the terms above the dimension vanish in the Chow ring, so drop them at each step
    M = variety(F)
    x = chern_character(F)
    out = x
    for _ in 2:n
      out = __truncate(M, out * x, dimension(M))
    end
    setfield!(new, :chern_character, out)
  end
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
  homog_chow = __has_chern_data(F) == __has_chern_data(G)
  !homog_chow && throw(ArgumentError("Dishomogeneous Chern data."))
  homog_weights = isdefined(F, :teleman_weights) == isdefined(G, :teleman_weights)
  !homog_weights && throw(ArgumentError("Dishomogeneous Teleman weights."))

  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, F.rank + G.rank)

  __has_chern_data(F) &&
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
  homog_chow = __has_chern_data(F) == __has_chern_data(G)
  !homog_chow && throw(ArgumentError("Dishomogeneous Chern data."))
  homog_weights = isdefined(F, :teleman_weights) == isdefined(G, :teleman_weights)
  !homog_weights && throw(ArgumentError("Dishomogeneous Teleman weights."))

  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, F.rank * G.rank)

  if __has_chern_data(F)
    # the terms above the dimension vanish in the Chow ring
    M = variety(F)
    setfield!(
      new,
      :chern_character,
      __truncate(M, chern_character(F) * chern_character(G), dimension(M)),
    )
  end
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
3*x2_1 + 3

julia> W = map(i -> exterior_power(F, i), 0:4);

julia> map(w -> (rank(w), chern_character(w)), W)
5-element Vector{Tuple{Int64, Singular.spoly{Singular.n_Q}}}:
 (1, 1)
 (3, 3*x2_1 + 3)
 (3, 6*x2_1 + 3)
 (1, 3*x2_1 + 1)
 (0, 0)
```
"""
function exterior_power(F::Bundle, k::Int)
  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, binomial(rank(F), k))

  __has_chern_data(F) &&
    setfield!(new, :chern_character, __simplify(__chern_characters_wedge(F, k)[end]))
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

The determinant is the line bundle with total Chern class ``1 + c_1(F)``, so it is
built from the first Chern class of `F` rather than as an exterior power.
"""
function det(F::Bundle)
  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, 1)

  if __has_chern_data(F)
    CH = F.parent.ring
    n = dimension(variety(F))
    classes = Dict{Int,Singular.spoly{Singular.n_Q}}(i => zero(CH) for i in 0:n)
    classes[0] = one(CH)
    n >= 1 && (classes[1] = chern_class(F, 1))
    setfield!(new, :chern_class, classes)
  end
  if isdefined(F, :teleman_weights)
    new_weights = Dict(
      hn_type => [sum(teleman_weights(F)[hn_type]; init=0)]
      for hn_type in keys(teleman_weights(F))
    )
    set_teleman_weights!(new, new_weights)
  end
  return new
end

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
3*x2_1 + 3

julia> W = map(i -> symmetric_power(F, i), 0:4);

julia> map(w -> (rank(w), chern_character(w)), W)
5-element Vector{Tuple{Int64, Singular.spoly{Singular.n_Q}}}:
 (1, 1)
 (3, 3*x2_1 + 3)
 (6, 12*x2_1 + 6)
 (10, 30*x2_1 + 10)
 (15, 60*x2_1 + 15)
```
"""
function symmetric_power(F::Bundle, k::Int)
  new = Bundle()
  setfield!(new, :parent, F.parent)
  setfield!(new, :rank, binomial(rank(F) + k - 1, rank(F) - 1))

  __has_chern_data(F) &&
    setfield!(
      new, :chern_character, __simplify(__chern_characters_symmetric(F, k)[end])
    )
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

"""
    __chern_characters_wedge(F::Bundle, k)

Compute the exterior powers of `F` up to degree `k`.
For internal use only.
"""
__chern_characters_wedge(F::Bundle, k) =
  __chern_characters_wedge(variety(F), chern_character(F), k)

"""
    __chern_characters_symmetric(F::Bundle, k)

Compute the symmetric powers of `F` up to degree `k`.
For internal use only.
"""
__chern_characters_symmetric(F::Bundle, k) =
  __chern_characters_symmetric(variety(F), chern_character(F), rank(F), k)

"""
    adams(F::Bundle, k)

Compute the Adams operation ``\\Phi^k`` on the Chern character of `F`.
For internal use only.
"""
adams(F::Bundle, k) = adams(variety(F), chern_character(F), k)

__chern_classes_from_character(F::Bundle) =
  __chern_classes_from_character(variety(F), chern_character(F))

__chern_character_from_classes(F::Bundle) =
  __chern_character_from_classes(variety(F), chern_classes(F), rank(F))

"""
    line_bundle(M::QuiverModuliSpace, eta::AbstractVector{Int}; unsafe::Bool=false, teleman::Bool=true)

Construct the descent of `L(eta)` on the quiver moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a quiver moduli space.
- `eta::AbstractVector{Int}`: a dimension vector.
- `unsafe::Bool`: Optional keyword argument to skip Chow ring computation checks.
  Default is `false`.
- `teleman::Bool`: Optional keyword argument to compute
  the Teleman weights of the line bundle. Default is `true`.

# Output

- the line bundle on `M` with Chern class `chern_class_line_bundle(M, eta)`.

# Examples

The ample generator of the Picard group of our favourite 6-fold:

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> chow_ring(M; chi=[-1, 1]);

julia> H = line_bundle(M, [3, -2]);

julia> chern_class(H)
-x2_1 + 1

julia> teleman_weights(H) ==
        Dict{HNType, Vector{Int64}}(
             [[2, 2], [0, 1]]         => [-10],
             [[2, 1], [0, 2]]         => [-20],
             [[1, 0], [1, 2], [0, 1]] => [-40],
             [[1, 0], [1, 3]]         => [-15],
             [[1, 0], [1, 1], [0, 2]] => [-35],
             [[1, 1], [1, 2]]         => [-5],
             [[2, 0], [0, 3]]         => [-30]
       )
true
```

The 7-subspace quiver moduli space:
```jldoctest
julia> Q = subspace_quiver(7); d = push!(ones(Int, 7), 2); M = QuiverModuliSpace(Q, d);

julia> K_dual = line_bundle(M, canonical_stability(Q, d); unsafe=true, teleman=false); K_dual.chern_class[1]
-2*x1_1 - 2*x2_1 - 2*x3_1 - 2*x4_1 - 2*x5_1 - 2*x6_1 + 7*x8_1
```
"""
function line_bundle(
  M::QuiverModuliSpace, eta::AbstractVector{Int}; unsafe::Bool=false, teleman::Bool=true
)
  eta' * M.d != 0 && throw(ArgumentError("$(eta) is not a linearization."))
  c1 = chern_class_line_bundle(M, eta; unsafe=unsafe)
  new = Bundle(M, 1, one(parent(c1)) + c1)
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
 2*x1_1 + 1
 3*x1_1 + 1
 4*x1_1 + 1
 5*x1_1 + 1
 6*x1_1 + 1

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
3*x2_1 + 1

julia> QuiverTools.teleman_weights(omega) ==
        Dict{HNType, Vector{Int64}}(
             [[2, 2], [0, 1]]         => [30],
             [[2, 1], [0, 2]]         => [60],
             [[1, 0], [1, 2], [0, 1]] => [120],
             [[1, 0], [1, 3]]         => [45],
             [[1, 0], [1, 1], [0, 2]] => [105],
             [[1, 1], [1, 2]]         => [15],
             [[2, 0], [0, 3]]         => [90]
       )
true
```
"""
function canonical_bundle(M::QuiverModuliSpace; teleman::Bool=true, unsafe::Bool=false)
  if !unsafe
    has_properly_semistable_representation(M.Q, M.d, M.theta, M.denom) &&
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
  end
  return line_bundle(M, -canonical_stability(M.Q, M.d); unsafe=unsafe, teleman=teleman)
end

"""
    universal_bundle(M:::QuiverModuliSpace, i::Int; teleman::Bool=true, unsafe::Bool=false)

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
 x1_1 + 1
 x2_1 + 1

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
 1//6*x1_2*x2_1 + 1//2*x2_1^2 + 1//12*x1_2*x2_2 - 1//36*x2_2^2 - 7//72*x2_1*x2_3 - 1//120*x2_2*x2_3 - 1//720*x2_3^2 - x1_2 + x2_1 - 1//2*x2_3 + 2
 2//3*x1_2*x2_1 + 1//2*x2_1^2 - 1//6*x1_2*x2_2 - 1//2*x2_1*x2_2 + 1//12*x2_2^2 + 1//24*x2_1*x2_3 + 1//180*x2_2*x2_3 + x2_1 - x2_2 + 3

julia> w = map(teleman_weights, [u1, u2]);

julia> w[1] ==
        Dict{HNType, Vector{Int64}}(
             [[2, 2], [0, 1]]         => [5, 5],
             [[2, 1], [0, 2]]         => [10, 10],
             [[1, 0], [1, 2], [0, 1]] => [25, 15],
             [[1, 0], [1, 3]]         => [10, 5],
             [[1, 0], [1, 1], [0, 2]] => [20, 15],
             [[1, 1], [1, 2]]         => [5, 0],
             [[2, 0], [0, 3]]         => [15, 15]
       )
true

julia> w[2] ==
        Dict{HNType, Vector{Int64}}(
             [[2, 2], [0, 1]]         => [5, 5, 0],
             [[2, 1], [0, 2]]         => [10, 5, 5],
             [[1, 0], [1, 2], [0, 1]] => [15, 15, 10],
             [[1, 0], [1, 3]]         => [5, 5, 5],
             [[1, 0], [1, 1], [0, 2]] => [15, 10, 10],
             [[1, 1], [1, 2]]         => [5, 0, 0],
             [[2, 0], [0, 3]]         => [10, 10, 10]
       )
true
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
5-element Vector{Singular.n_Q}:
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
2*x3_1 + 1

julia> degree(F)
56
```

The 7-subspace quiver:
```jldoctest
julia> Q = subspace_quiver(7); d = push!(ones(Int, 7), 2); M = QuiverModuliSpace(Q, d);

julia> K_dual = line_bundle(M, canonical_stability(Q, d); unsafe=true);

julia> degree(K_dual; unsafe=true)
154
```
"""
function degree(F::Bundle; unsafe::Bool=false)
  M = variety(F)
  unsafe ? (n = 1 - euler_form(M.Q, M.d, M.d)) : (n = dimension(M))

  if rank(F) == 1
    out = chern_class(F)
  else
    out = chern_class(det(F))
  end

  out = __simplify(out)

  if n >= 2
    class_det = deepcopy(out)
    for _ in 1:(n - 1)
      # the multiplication here takes most of runtime
      # what about out *= class_det?
      Singular.mul!(out, out, class_det)
      out = __simplify(out)
    end
  end

  return integral(M, out; unsafe=unsafe)
end

function __format_chern_monomial(partition::AbstractVector{Int})
  isempty(partition) && return "1"
  pieces = map(unique(partition)) do k
    m = count(==(k), partition)
    m == 1 ? "c_$k" : "c_$k^$m"
  end
  return join(pieces, " ")
end

function __format_universal_monomial(
  generator_labels::AbstractVector{Tuple{Int,Int}}, exponents::AbstractVector{Int}
)
  # the empty monomial is the fundamental class
  all(iszero, exponents) && return "1"
  factors = String[]
  for g in eachindex(exponents)
    exponents[g] == 0 && continue
    vertex, chern_degree = generator_labels[g]
    symbol = "c_$(chern_degree)(U_$(vertex))"
    push!(factors, exponents[g] == 1 ? symbol : "$(symbol)^$(exponents[g])")
  end
  return join(factors, " ")
end

"""
    __weighted_monomials(weights::AbstractVector{Int}, n::Int)

Enumerate the exponent vectors of all monomials of weighted degree `n` in
generators whose individual degrees are given by `weights`. Used to list every
monomial of a fixed codimension in a graded ring. For internal use only.
"""
function __weighted_monomials(degrees::AbstractVector{Int}, total_degree::Int)
  monomials = Vector{Int}[]
  exponents = zeros(Int, length(degrees))
  __extend_weighted_monomials!(monomials, exponents, degrees, 1, total_degree)
  return monomials
end

# assign an exponent to generator `gen`, then recurse on the remaining
# generators, spending `remaining_degree` of the total budget across them.
# This is a top-level function rather than a closure inside `__weighted_monomials`,
# because a recursive closure is stored in a `Core.Box` and every call to it dispatches
# at run time.
function __extend_weighted_monomials!(
  monomials::Vector{Vector{Int}},
  exponents::Vector{Int},
  degrees::AbstractVector{Int},
  gen::Int,
  remaining_degree::Int,
)
  if gen > length(degrees)
    # all exponents assigned; keep this monomial iff the budget is exactly spent
    remaining_degree == 0 && push!(monomials, copy(exponents))
    return nothing
  end
  for exponent in 0:(remaining_degree ÷ degrees[gen])
    exponents[gen] = exponent
    __extend_weighted_monomials!(
      monomials, exponents, degrees, gen + 1, remaining_degree - exponent * degrees[gen]
    )
  end
  exponents[gen] = 0  # reset before unwinding to the caller
  return nothing
end

"""
    __degree1_vector(nf, L::Int)

Read the coordinate vector (over the degree-1 variables of the Chow ring) of a
homogeneous codimension-1 class from its Gröbner normal form `nf`. For internal
use only.
"""
function __degree1_vector(normal_form, n_variables::Int)
  coordinates = zeros(Rational{BigInt}, n_variables)
  # a linear form is a sum of terms `coefficient * x_k`; record each coefficient
  # at the index of the single variable appearing in its term
  for (coefficient, exponent_vector) in
      zip(Singular.coefficients(normal_form), Singular.exponent_vectors(normal_form))
    var_index = findfirst(==(1), exponent_vector)
    var_index === nothing && continue  # skip any constant/non-linear term
    coordinates[var_index] = Rational{BigInt}(
      BigInt(Singular.numerator(coefficient)), BigInt(Singular.denominator(coefficient))
    )
  end
  return coordinates
end

"""
    __universal_chern_generators(M::QuiverModuliSpace)

Compute the universal Chern classes ``c_j(U_i)`` that minimally generate the Chow
ring `chow_ring(M)`.
Return triples `(vertex, chern_degree, var_index)` where `var_index` is the
index of the generator in `gens(CH)`.

The first Chern classes ``c_1(U_i)`` satisfy linear relations (from the choice of
linearization), so only a maximal linearly independent subset of them is kept;
they are processed in vertex order, so that ``c_1(U_1)`` is the first.
All higher Chern classes are always kept. For internal use only.
"""
function __universal_chern_generators(M::QuiverModuliSpace)
  generators = gens(chow_ring(M))
  n_variables = length(generators)
  # row-echelon basis of the kept first Chern classes, used to test the linear
  # independence of each new c_1(U_i) against those already kept
  echelon_rows = Vector{Vector{Rational{BigInt}}}()
  kept = Tuple{Int,Int,Int}[]
  var_index = 0

  # Iterate over each generator of the Chow ring,
  # compute its coordinates,
  # perform row reduction with the previously kept first Chern classes
  # and if it is linearly independent add it to the echelon basis.
  for vertex in 1:n_vertices(M.Q), chern_degree in 1:M.d[vertex]
    var_index += 1
    # higher Chern classes are independent generators of the Chow ring
    if chern_degree > 1
      push!(kept, (vertex, chern_degree, var_index))
      continue
    end

    coords = __degree1_vector(__simplify(generators[var_index]), n_variables)

    for row in echelon_rows
      pivot = findfirst(!iszero, row)
      iszero(coords[pivot]) || (coords .-= (coords[pivot] / row[pivot]) .* row)
    end
    # if linearly independent, we keep this generator.
    if any(!iszero, coords)
      push!(echelon_rows, coords)
      sort!(echelon_rows; by=row -> findfirst(!iszero, row))
      push!(kept, (vertex, chern_degree, var_index))
    end
  end
  return kept
end

"""
    chern_numbers(F::Bundle)

Compute the Chern numbers of a bundle `F`, i.e. all top intersection products
``\\int_M c_{i_1}(F) \\cdots c_{i_k}(F)``, indexed by partitions `[i_1, ..., i_k]` of `dimension(M)`.

Returns a `Dict{String,BigInt}`, since the Chern numbers of high-dimensional moduli
spaces do not fit in an `Int`.
"""
function chern_numbers(F::Bundle; unsafe::Bool=false)
  M = variety(F)
  CH = chow_ring(M; unsafe=unsafe)
  classes = chern_classes(F)

  # the Chern numbers of an element in the Chow ring
  # are the top intersection product c_{i_1} ... c_{i_k}
  # per every partition (i_1, ..., i_k) of the dimension of the underlying variety.
  return Dict{String,BigInt}(
    __format_chern_monomial(partition) => BigInt(
      Singular.numerator(
        integral(M, prod(classes[k] for k in partition; init=one(CH)); unsafe=unsafe)
      ),
    )
    for partition in partitions(dimension(M))
  )
end

"""
    chern_numbers(M::QuiverModuliSpace; unsafe::Bool=false, universal::Bool=false)

Compute the Chern numbers of the tangent bundle of the quiver moduli space `M`,
i.e. all top intersection products
``\\int_M c_{i_1}(T_M) \\cdots c_{i_k}(T_M)``,
indexed by partitions `[i_1, ..., i_k]` of `dimension(M)`.

Returns a `Dict{String,BigInt}` whose keys describe each monomial in the Chern
classes and whose values are the corresponding Chern numbers.

If `universal` is set to `true`, the result additionally contains all top
intersection products of the Chern classes ``c_j(U_i)`` of the universal
bundles, which generate the Chow ring (cf [[Theorem 14, doi:10.1007/s00209-015-1409-8](https://doi.org/10.1007/s00209-015-1409-8)].
These keys are of the form `c_j(U_i)`,
e.g. `c_2(U_1)` or `c_2(U_1) c_3(U_2)`.

The first Chern classes ``c_1(U_i)`` are linearly dependent (through the chosen
linearization), so only a minimal generating set is used: a maximal linearly
independent subset of them is kept, favoring ``c_1(U_1)``. All higher Chern
classes are kept.

# Example

The Chern numbers of `P^3`, realized as the moduli space of the
generalized Kronecker quiver with 4 arrows and dimension vector `(1, 1)`:

```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(4), [1, 1]); chow_ring(M);

julia> cn = chern_numbers(M);

julia> cn["c_1^3"], cn["c_2 c_1"], cn["c_3"]
(64, 24, 4)

julia> cn = chern_numbers(M; universal=true);

julia> @assert cn == Dict(
             "c_1^3" => 64,
             "c_2 c_1" => 24,
             "c_3" => 4,
             "c_1(U_1)^3" => -1,
             )

```


The Chern numbers of our favourite 6-fold:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> cn = chern_numbers(M);

julia> @assert cn == Dict(
               "c_1^6" => 41553,
               "c_2 c_1^4" => 21141,
               "c_2^2 c_1^2" => 10827,
               "c_2^3" => 5569,
               "c_3 c_1^3" => 6561,
               "c_3 c_2 c_1" => 3357,
               "c_3^2" => 1041,
               "c_4 c_1^2" => 1269,
               "c_4 c_2" => 643,
               "c_5 c_1" => 153,
               "c_6" => 13,
              )

```

We compute the Chern numbers of the top intersection products of the Chern classes
of the universal bundles as well, replicating the result of [[Proposition 3.8, arXiv:2412.15390](https://arxiv.org/pdf/2412.15390)].
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> cn = chern_numbers(M; universal=true);

julia> universal = Dict(k => v for (k, v) in cn if occursin("(U_", k));

julia> @assert universal == Dict(
             "c_1(U_1) c_2(U_1) c_3(U_2)" => 2,
             "c_1(U_1) c_2(U_2) c_3(U_2)" => 3,
             "c_1(U_1)^2 c_2(U_1) c_2(U_2)" => 9,
             "c_1(U_1)^2 c_2(U_1)^2" => 6,
             "c_1(U_1)^2 c_2(U_2)^2" => 14,
             "c_1(U_1)^3 c_3(U_2)" => 5,
             "c_1(U_1)^4 c_2(U_1)" => 18,
             "c_1(U_1)^4 c_2(U_2)" => 27,
             "c_1(U_1)^6" => 57,
             "c_2(U_1) c_2(U_2)^2" => 5,
             "c_2(U_1)^2 c_2(U_2)" => 3,
             "c_2(U_1)^3" => 2,
             "c_2(U_2)^3" => 9,
             "c_3(U_2)^2" => 1,
            )

```



The `unsafe` keyword argument speeds up computations,
skipping checks for ample stability and existence
of properly semistable representations.
We showcase its use in the Chern numbers
of the 7-subspace quiver moduli space:

```jldoctest
julia> Q = subspace_quiver(7); d = push!(ones(Int, 7), 2);

julia> M = QuiverModuliSpace(Q, d);

julia> cn = chern_numbers(M; unsafe=true);

julia> @assert cn == Dict(
         "c_1^4" => 154,
         "c_2 c_1^2" => 112,
         "c_2^2" => 136,
         "c_3 c_1" => 56,
         "c_4" => 38,
       )

```
"""
function chern_numbers(M::QuiverModuliSpace; unsafe::Bool=false, universal::Bool=false)
  CH = chow_ring(M; unsafe=unsafe)
  T = tangent_bundle(M; unsafe=unsafe)

  out = chern_numbers(T; unsafe=unsafe)
  !universal && return out

  # the `universal` keyword for now is a convenient way to compute
  # top intersection products of Chern classes of the universal bundles specifically.
  generators = gens(CH)

  kept_generators_indices = __universal_chern_generators(M)
  # the actual indices of the variables we need
  generator_indices = Int[var_index for (_, _, var_index) in kept_generators_indices]
  # to label the output, we keep track of vertices and degrees.
  generator_labels = Tuple{Int,Int}[
    (vertex, chern_degree) for (vertex, chern_degree, _) in kept_generators_indices
  ]
  # to enumerate all the top degree monomials.
  generator_degrees = Int[chern_degree for (_, chern_degree) in generator_labels]

  for exponents in __weighted_monomials(generator_degrees, dimension(M))
    # the monomial is a single term in the variables of the Chow ring, so build it directly
    full_exponents = zeros(Int, length(generators))
    full_exponents[generator_indices] = exponents
    context = Singular.MPolyBuildCtx(CH)
    Singular.push_term!(context, one(Singular.base_ring(CH)), full_exponents)
    monomial = Singular.finish(context)

    out[__format_universal_monomial(generator_labels, exponents)] = BigInt(
      Singular.numerator(integral(M, monomial; unsafe=unsafe))
    )
  end

  return out
end
