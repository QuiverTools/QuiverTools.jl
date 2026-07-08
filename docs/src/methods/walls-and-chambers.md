# Walls and Chambers decomposition

QuiverTools implements the walls-and-chambers decomposition
of the GIT problem of quiver moduli using the polyhedral geometry
interface of `Oscar.jl`.

This functionality lives in a package extension and is only available once Oscar
is loaded. Load it with `import Oscar` (rather than `using Oscar`, which would pull
Oscar's exports into scope and clash with QuiverTools names such as `index` and
`todd_class`):

```julia
using QuiverTools
import Oscar
```

Without Oscar loaded, the functions below raise an error explaining this.

```@docs
is_special_subdimension_vector
all_special_subdimension_vectors
sst
vgit_walls
wall_system
vgit_chambers
vgit_fan
git_equivalent
```
