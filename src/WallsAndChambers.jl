# Walls-and-chambers and VGIT functionality relies on Oscar's polyhedral geometry
# and therefore lives in the Oscar package extension (ext/QuiverToolsOscarExt.jl).
# The functions below are only method stubs: the extension adds the real methods
# once Oscar is loaded. Until then, calling any of them raises a clear error.
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

```jldoctests
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

```jldoctests
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

```jldoctests
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

```jldoctests
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

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-4"); d = [1, 1, 1, 1];

julia> map(Oscar.rays, vgit_chambers(Q, d; verbose=false))
3-element Vector{Oscar.SubObjectIterator{Oscar.RayVector{Nemo.QQFieldElem}}}:
 [[1, -1, 0, 0], [0, 0, 1, -1], [1, 0, 0, -1]]
 [[0, 1, -1, 0], [0, 0, 1, -1], [1, 0, 0, -1]]
 [[1, -1, 0, 0], [1, 0, 0, -1], [0, 1, -1, 0]]
```

The following example has three inner walls as well, but one is equal to its wall system
hyperplane and the two others are not.

```jldoctests
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

```jldoctests
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

Requires Oscar: run `import Oscar` to enable this function.

# Example

```jldoctests
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
```jldoctests
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
function git_equivalent end

# Load Oscar on demand. This pulls in QuiverToolsOscarExt, which defines the
# concrete methods that shadow the catch-all stubs below. We `import` rather than
# `using` Oscar so its exports do not clash with QuiverTools names such as `index`
# and `todd_class`; the extension triggers on either.
function _load_oscar()
  Base.get_extension(@__MODULE__, :QuiverToolsOscarExt) === nothing || return nothing
  @eval Main import Oscar
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
  :vgit_walls, :wall_system, :vgit_chambers, :vgit_fan, :git_equivalent,
)
  @eval @oscar_stub $f
end
