# Walls and Chambers decomposition

QuiverTools implements the walls-and-chambers decomposition
of the GIT problem of quiver moduli.
The Boolean equivalence test [`git_equivalent`](@ref) is implemented combinatorially
and has no optional dependencies.
Constructing the polyhedral cones and fans uses the geometry interface of `Oscar.jl`.

The polyhedral functionality lives in a *package extension*
that depends on the `Oscar` algebra system.

To use it, one must

- Ensure `Oscar` is installed in the current environment, by running `using Pkg; Pkg.add("Oscar")`; and
- Load the Oscar package, possibly without its interface, by running `import Oscar`.

If `import Oscar` is not run, `QuiverTools` will run it automatically,
but this may trigger recompilation and produce unwanted output.
For non-interactive workflows, it is best to run it at the beginning of the session.

Note that running `using Oscar` instead of `import Oscar` will expose
the Oscar interface, and this might clash with functions provided by QuiverTools.
The interface of either package can be accessed even if it is not exposed,
by prepending the package name to the bindings. 
For instance, if one only runs `import Oscar`,
the function `dim()` of Oscar can be used
by running `Oscar.dim()`.


```@docs
is_special_subdimension_vector
all_special_subdimension_vectors
sst
vgit_walls
wall_system
vgit_chambers
vgit_fan
git_equivalent
all_stability_parameters
```
