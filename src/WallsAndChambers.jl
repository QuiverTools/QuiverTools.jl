# The combinatorial VGIT layer is independent of Oscar. It records the equations and
# inequalities defining semistable cones and walls, and provides exact point and line
# segment predicates. The Oscar extension consumes the same data to construct the
# corresponding polyhedra.

# Equations and inequalities defining sst(d). The inequalities are Schofield's general
# subdimension vectors; zero and d impose no inequalities and are omitted.
function __semistable_cone_data(Q::Quiver, d::AbstractVector{Int})
  return __semistable_cone_data(Q, Vector{Int}(d))
end

@memoize Dict function __semistable_cone_data(Q::Quiver, d::Vector{Int})
  inequalities = filter(
    e -> any(!=(0), e) && e != d,
    all_general_subdimension_vectors(Q, d),
  )
  return (equations=[d], inequalities=inequalities)
end

# Equations and inequalities defining W_e = sst(e) ∩ sst(d-e) ∩ sst(d).
function __vgit_wall_data(
  Q::Quiver,
  d::AbstractVector{Int},
  e::AbstractVector{Int},
)
  return __vgit_wall_data(Q, Vector{Int}(d), Vector{Int}(e))
end

@memoize Dict function __vgit_wall_data(Q::Quiver, d::Vector{Int}, e::Vector{Int})
  data = __semistable_cone_data.(Ref(Q), (e, d - e, d))
  equations = unique!(vcat((datum.equations for datum in data)...))
  inequalities = unique!(vcat((datum.inequalities for datum in data)...))
  return (equations=equations, inequalities=inequalities)
end

# Evaluate a scalar product without losing exactness when the parameter is rational.
function __exact_dot(x::AbstractVector{Int}, y::AbstractVector)
  length(x) == length(y) || throw(DimensionMismatch("vectors must have equal lengths"))
  value = big(0) // big(1)
  for i in eachindex(x, y)
    value += big(x[i]) * y[i]
  end
  return value
end

function __in_rational_cone(data, theta::AbstractVector)
  return all(normal -> iszero(__exact_dot(normal, theta)), data.equations) &&
         all(normal -> __exact_dot(normal, theta) <= 0, data.inequalities)
end

"""Return whether `theta` belongs to the semistable cone `sst(d)`."""
function __in_semistable_cone(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector,
)
  return __in_rational_cone(__semistable_cone_data(Q, d), theta)
end

"""Return whether `theta` belongs to the VGIT wall `W_e`."""
function __in_vgit_wall(
  Q::Quiver,
  d::AbstractVector{Int},
  e::AbstractVector{Int},
  theta::AbstractVector,
)
  return __in_rational_cone(__vgit_wall_data(Q, d, e), theta)
end

# Intersect the segment theta1--theta2 with a rational polyhedral cone. The result is
# the exact closed interval of parameters t in [0,1] for which
# (1-t)theta1 + t theta2 belongs to the cone, or nothing when it is empty.
function __segment_cone_intersection(data, theta1::AbstractVector, theta2::AbstractVector)
  length(theta1) == length(theta2) ||
    throw(DimensionMismatch("stability parameters must have equal lengths"))
  lower = big(0) // big(1)
  upper = big(1) // big(1)

  for normal in data.equations
    initial = __exact_dot(normal, theta1)
    delta = __exact_dot(normal, theta2) - initial
    if iszero(delta)
      iszero(initial) || return nothing
    else
      crossing = -initial / delta
      lower = max(lower, crossing)
      upper = min(upper, crossing)
      lower <= upper || return nothing
    end
  end

  for normal in data.inequalities
    initial = __exact_dot(normal, theta1)
    delta = __exact_dot(normal, theta2) - initial
    if iszero(delta)
      initial <= 0 || return nothing
    elseif delta > 0
      upper = min(upper, -initial / delta)
    else
      lower = max(lower, -initial / delta)
    end
    lower <= upper || return nothing
  end

  return (lower, upper)
end

"""Return whether `theta` lies in a VGIT chamber rather than on a wall."""
function __is_vgit_chamber_parameter(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector,
)
  __in_semistable_cone(Q, d, theta) || return false
  return all(all_subdimension_vectors(d; nonzero=true, strict=true)) do e
    !__in_vgit_wall(Q, d, e, theta)
  end
end

# The target may meet a wall at the endpoint of the segment, but the segment must not
# meet any wall earlier. Using complete cone data also handles smaller walls when the
# whole segment lies in their defining hyperplane.
function __in_closure_of_vgit_chamber(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector,
  thetabar::AbstractVector,
)
  __is_vgit_chamber_parameter(Q, d, theta) || return false
  __in_semistable_cone(Q, d, thetabar) || return false
  endpoint = big(1) // big(1)
  for e in all_subdimension_vectors(d; nonzero=true, strict=true)
    intersection = __segment_cone_intersection(__vgit_wall_data(Q, d, e), theta, thetabar)
    intersection === nothing && continue
    intersection == (endpoint, endpoint) || return false
  end
  return true
end

# Polyhedral wall and chamber objects rely on Oscar and therefore live in the Oscar
# package extension (ext/QuiverToolsOscarExt.jl). The declarations below are method
# stubs which load that extension on demand.
#
# To enable them, load Oscar alongside QuiverTools. Use `import Oscar` rather than
# `using Oscar`: it activates the extension without bringing Oscar's exports into
# scope, which would clash with QuiverTools names such as `index` and `todd_class`.
#
#     using QuiverTools
#     import Oscar

"""
    is_special_subdimension_vector(Q::Quiver, e::AbstractVector{Int}, d::AbstractVector{Int})

Compute whether `e` is a special subdimension vector of `d` for `Q`.

Special subdimension vectors are defined in
[[Definition 6.1, MR5007902](https://mathscinet.ams.org/mathscinet-getitem?mr=5007902)].

Requires Oscar: run `import Oscar` to enable this function.

# Example

```jldoctest
julia> Q = Quiver("1-2,2-3,3-4,1-4"); d = [1, 1, 1, 1];

julia> is_special_subdimension_vector(Q, [1, 1, 1, 1], d)
false

julia> is_special_subdimension_vector(Q, [1, 0, 0, 1], d)
true
```
"""
function is_special_subdimension_vector end

"""
    all_special_subdimension_vectors(Q::Quiver, d::AbstractVector{Int})

Compute all the special subdimension vectors of `d` for the quiver `Q`.

Requires Oscar: run `import Oscar` to enable this function.

# Example

```jldoctest
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> length(all_special_subdimension_vectors(Q, d))
7
```
"""
function all_special_subdimension_vectors end

"""
    sst(Q, e)

Compute the semistable cone `sst(e)` for a given quiver `Q` and a vector `e`.
This is the cone of all stability parameters for which semistable representations
of the quiver `Q` with dimension vector `e` exist.

Requires Oscar: run `import Oscar` to enable this function.

# Example

```jldoctest
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> collect(Oscar.rays(sst(Q, d)))
3-element Vector{Oscar.RayVector{Nemo.QQFieldElem}}:
 [0, 0, 1, -1]
 [0, 1, -1, 0]
 [1, -1, 0, 0]
```
"""
function sst end

"""
    vgit_walls(Q, d; inner=false, top_dimension=true)

Compute all the walls `W_e` of the quiver `Q` with dimension vector `d`.
Defaults to only computing the top-dimensional walls, pass
the `top_dimension=false` keyword to compute all the `W_e` of the VGIT problem.
Defaults to computing all walls, pass the `inner=true` keyword
to not include the outer walls.

Requires Oscar: run `import Oscar` to enable this function.

# Example

```jldoctest
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> walls = vgit_walls(Q, d; inner=false, top_dimension=false);

julia> map(Oscar.dim, walls)
7-element Vector{Int64}:
 2
 2
 2
 2
 1
 2
 2
```
"""
function vgit_walls end

"""
    wall_system(Q, d; inner=false, as_cones=true)

Compute the wall system for the quiver `Q` with dimension vector `d`.
This is the set of hyperplanes `H_e` for which there exists at least one
vgit wall `W_e` that lays on `H_e`.
By default it returns the walls intersected with the semistable cone `sst(d)`,
pass the keyword `as_cones=false` to get the full hyperplanes.

Requires Oscar: run `import Oscar` to enable this function.
"""
function wall_system end

"""
    vgit_chambers(Q, d; verbose=false)

Compute all VGIT chambers for the quiver `Q` with dimension vector `d`.

VGIT chambers are the top-dimensional equivalence classes in the VGIT problem;
if stable representations exist,
then the VGIT chambers have dimension equal to `Oscar.dim(sst(Q, d))`
and are defined by the walls `W_e` of codimension 1.

If there exit no stable representations,
`Oscar.dim(sst(Q, d))` is strictly smaller than `length(d) - 1`.

Requires Oscar: run `import Oscar` to enable this function.

# Example

The following example has three inner walls `W_e`, and all of them are strictly smaller
than the corresponding `H_e \\cap sst(d)` in the wall system.

```jldoctest
julia> Q = Quiver("1-2,2-3,3-4,1-4"); d = [1, 1, 1, 1];

julia> map(Oscar.rays, vgit_chambers(Q, d; verbose=false))
3-element Vector{Oscar.SubObjectIterator{Oscar.RayVector{Nemo.QQFieldElem}}}:
 [[1, -1, 0, 0], [0, 0, 1, -1], [1, 0, 0, -1]]
 [[0, 1, -1, 0], [0, 0, 1, -1], [1, 0, 0, -1]]
 [[1, -1, 0, 0], [1, 0, 0, -1], [0, 1, -1, 0]]
```

The following example has three inner walls as well, but one is equal to its wall system
hyperplane and the two others are not.

```jldoctest
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> map(Oscar.rays, vgit_chambers(Q, d; verbose=false))
4-element Vector{Oscar.SubObjectIterator{Oscar.RayVector{Nemo.QQFieldElem}}}:
 [[1, -1, 0, 0], [1, 0, 0, -1], [1, 0, -1, 0]]
 [[1, -1, 0, 0], [0, 0, 1, -1], [1, 0, 0, -1]]
 [[0, 1, -1, 0], [1, 0, 0, -1], [1, 0, -1, 0]]
 [[0, 1, -1, 0], [0, 0, 1, -1], [1, 0, 0, -1]]
```
"""
function vgit_chambers end

"""
    vgit_fan(Q, d; verbose=false)

Compute the VGIT fan for the quiver `Q` with dimension vector `d`.

Requires Oscar: run `import Oscar` to enable this function.

# Example

```jldoctest
julia> Q = three_vertex_quiver(2, 3, 4); d = [2, 3, 4];

julia> F = vgit_fan(Q, d); Oscar.rays(F)
9-element Oscar.SubObjectIterator{Oscar.RayVector{Nemo.QQFieldElem}}:
 [1, -2//3, 0]
 [1, -1//2, -1//8]
 [1, -2//5, -1//5]
 [1, 0, -1//2]
 [1, 2//9, -2//3]
 [1, 2//5, -4//5]
 [1, 2//3, -1]
 [1, 2, -2]
 [0, 1, -3//4]
```
"""
function vgit_fan end

"""
    git_equivalent(Q, d, theta1, theta2)

Check if the two stability parameters `theta1` and `theta2` are equivalent.

By [[Corollary 4.4, MR5007902](https://mathscinet.ams.org/mathscinet-getitem?mr=5007902)],
this is equivalent to their convex hull
either lying in a wall or not intersecting any of them.

This computation uses the combinatorial equations and inequalities defining the VGIT
walls and does not require Oscar.

# Example

```jldoctest
julia> Q = three_vertex_quiver(2, 3, 4); d = [1, 2, 2];

julia> F = vgit_fan(Q, d); Oscar.rays(F)
4-element Oscar.SubObjectIterator{Oscar.RayVector{Nemo.QQFieldElem}}:
 [1, -1//2, 0]
 [1, 0, -1//2]
 [1, 1//2, -1]
 [0, 1, -1]

julia> theta1 = [2, -1//2, -1//2];

julia> theta2 = [2, 1//2, -3//2];

julia> theta3 = [1, 3//2, -2];

julia> git_equivalent(Q, d, theta1, theta2)
false

julia> git_equivalent(Q, d, theta1, theta3)
false

julia> git_equivalent(Q, d, theta2, theta3)
false
```

An example where one stability parameter lies on a wall:
```jldoctest
julia>  Q = Quiver("1---------2,1-3,2---3"); d = [1, 2, 3];

julia> x = [5, -1, -1]; y = [3, 0, -1];

julia> git_equivalent(Q, d, x, y)
false

julia> # indeed,

julia> W = vgit_walls(Q, d; top_dimension=false);

julia> any(x in w for w in W)
false

julia> any(y in w for w in W)
true
```
"""
function git_equivalent(
  Q::Quiver,
  d::AbstractVector{Int},
  theta1::AbstractVector,
  theta2::AbstractVector,
)
  theta1 == theta2 && return true
  for e in all_subdimension_vectors(d; nonzero=true, strict=true)
    data = __vgit_wall_data(Q, d, e)
    __segment_cone_intersection(data, theta1, theta2) === nothing && continue
    __in_rational_cone(data, theta1) && __in_rational_cone(data, theta2) && continue
    return false
  end
  return true
end

"""
    all_stability_parameters(Q::Quiver, d::AbstractVector{Int}; generic::Bool=false)

Compute a list of stability parameters, one for each equivalence class of
stability parameters for the quiver `Q` and dimension vector `d`.

For now this excludes the semisimple condition,
which corresponds to the stability parameter `zeros(Int, length(d))`.

The `generic` keyword argument only returns stability parameters
from the top-dimensional chambers of the VGIT fan.

Requires Oscar: run `import Oscar` to enable this function.

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> all_stability_parameters(Q, d)
17-element Vector{Vector{Int64}}:
 [1, -1, 0, 0]
 [1, 0, -1, 0]
 [1, 0, 0, -1]
 [0, 0, 1, -1]
 [0, 1, -1, 0]
 [2, -1, -1, 0]
 [2, -1, 0, -1]
 [2, 0, -1, -1]
 [1, 0, 1, -2]
 [1, -1, 1, -1]
 [1, 1, -2, 0]
 [1, 1, -1, -1]
 [0, 1, 0, -1]
 [3, -1, -1, -1]
 [2, -1, 1, -2]
 [2, 1, -2, -1]
 [1, 1, 0, -2]
```

The `generic` keyword argument only returns the stability parameters
from the top-dimensional chambers of the VGIT fan, i.e., the ones
for which semistability and stability are equivalent:

```jldoctest
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> all_stability_parameters(Q, d; generic=true)
4-element Vector{Vector{Int64}}:
 [3, -1, -1, -1]
 [2, -1, 1, -2]
 [2, 1, -2, -1]
 [1, 1, 0, -2]
```

This method behaves well with respect to the trivial case of our favourite quiver, of course:

```jldoctests
julia> Q = kronecker_quiver(3); d = [2, 3];

julia> all_stability_parameters(Q, d)
1-element Vector{Vector{Int64}}:
 [3, -2]
```
"""
function all_stability_parameters end

# Load Oscar on demand. This pulls in QuiverToolsOscarExt, which defines the
# concrete methods that shadow the catch-all stubs below. We `import` rather than
# `using` Oscar so its exports do not clash with QuiverTools names such as `index`
# and `todd_class`; the extension triggers on either.
function _load_oscar()
  Base.get_extension(@__MODULE__, :QuiverToolsOscarExt) === nothing || return nothing
  redirect_stdout(devnull) do
    @eval Main import Oscar
  end
  return nothing
end

"""
    @oscar_stub f

Mark `f` as an Oscar-backed entry point. Until Oscar is loaded, calling `f` loads
it on demand and re-dispatches to the concrete method the extension defines. Once
loaded, that method (being more specific than this varargs catch-all) is hit
directly, so the trigger only fires once. If Oscar is already loaded and still
nothing matches, this errors instead of looping.
"""
macro oscar_stub(f)
  quote
    function $(esc(f))(args...; kwargs...)
      if Base.get_extension(@__MODULE__, :QuiverToolsOscarExt) !== nothing
        throw(MethodError($(esc(f)), args))
      end
      _load_oscar()
      return Base.invokelatest($(esc(f)), args...; kwargs...)
    end
  end
end

# When Oscar is not loaded, every entry point above resolves to its catch-all,
# which loads Oscar and re-dispatches. The extension adds concrete methods that
# take precedence once Oscar has been loaded.
for f in (
  :is_special_subdimension_vector, :all_special_subdimension_vectors, :sst,
  :vgit_walls, :wall_system, :vgit_chambers, :vgit_fan, :all_stability_parameters,
)
  @eval @oscar_stub $f
end
