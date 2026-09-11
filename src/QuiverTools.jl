module QuiverTools

using Pkg

using Memoization: Memoization
using IterTools: IterTools
using LinearAlgebraX: LinearAlgebraX
using Combinatorics
using StaticArrays

using Singular: Singular

using PrecompileTools: PrecompileTools

import Base:
  show, ==, hash, convert, getindex, setindex!, length, iterate, keys, haskey, *, +, -, ^
import Memoization: @memoize, empty_all_caches!, empty_cache!
import IterTools: subsets
import LinearAlgebraX: rankx
import Combinatorics: combinations, with_replacement_combinations, partitions, permutations
import Singular: polynomial_ring, degree, coeff, constant_coefficient,
  preimage, quotient_ideal, fraction_field, std, gens, base_ring

# optional dependencies

# try
#   using SchubertPolynomials: SchubertPolynomials
#   import SchubertPolynomials: xy_ring, schub_poly
# catch # will not run as long as SchubertPolynomials is installed or in the [deps]
#   @warn "SchubertPolynomials not found. Chow ring default functionnality may be slower.
# To solve this, install SchubertPolynomials.jl by running\n\n
# using Pkg; Pkg.rm(\"QuiverTools\"); Pkg.add(url=\"https://github.com/pseudoeffective/SchubertPolynomials.jl\", rev=\"9f1d979b20860ccba474775e02cd15663f31e304\")\n\n
# Then force recompilation of QuiverTools.jl by running\n\n
# Base.compilecache(Base.identify_package("QuiverTools"))\n\n"
# end

# Types
export Quiver, HNType, LunaType, QuiverModuli, QuiverModuliSpace, QuiverModuliStack, Bundle

# Quivers
export n_vertices,
  n_arrows, arrows, indegree, outdegree, is_acyclic, is_connected, is_sink, is_source,
  strongly_connected_components, underlying_graph, first_hochschild_cohomology

# Constructors
export kronecker_quiver, loop_quiver, jordan_quiver, subspace_quiver, star_quiver,
  generalized_subspace_quiver, thickened_subspace_quiver, three_vertex_quiver,
  cyclic_quiver,
  bipartite_quiver, opposite_quiver, double_quiver, disjoint_union, dynkin_quiver,
  extended_dynkin_quiver
export kronecker_moduli, subspace_quiver_moduli

# Stability
export canonical_stability, is_coprime, slope
export all_hn_types,
  is_hn_type, has_semistables, has_stables, has_properly_semistables,
  codimension_hn_stratum, is_amply_stable
export is_general_subdimension_vector, all_general_subdimension_vectors

# Representation theory
export euler_form, euler_matrix, is_root, is_schur_root, is_real_root, is_imaginary_root,
  is_isotropic_root,
  general_ext, general_hom, canonical_decomposition, in_fundamental_domain,
  bocklandt_reduction, is_coregular, is_cofree

# Moduli
export all_luna_types, is_luna_type, dimension_of_luna_stratum
export is_nonempty, codimension_unstable_locus, codimension_singular_locus, dimension,
  is_smooth,
  is_projective, is_strongly_amply_stable, semistable_equals_stable, semisimple_moduli_space

# Hodge
export hodge_diamond, hodge_polynomial, picard_rank, index, betti_numbers

# Chow
export chow_ring, motive, index, betti_numbers, poincare_polynomial, is_smooth,
  is_projective,
  semisimple_moduli_space, point_class, todd_class, chern_class_line_bundle,
  chern_character_line_bundle, total_chern_class_universal,
  integral, euler_characteristic

# Teleman
export teleman_bounds, weights_hn_type, weights_universal_bundle, weights_canonical_bundle,
  all_weights_endomorphisms_universal_bundle, weights_endomorphisms_universal_bundles,
  does_rigidity_inequality_hold, does_teleman_inequality_hold, set_teleman_weights!

# Bundles
export chern_character, chern_class, chern_classes, dual, exterior_power, symmetric_power,
  det, line_bundle, canonical_bundle, universal_bundle, tangent_bundle, chern_numbers,
  degree, rank, teleman_weights, structure_sheaf

# Walls and Chambers
export rays, is_special_subdimension_vector, all_special_subdimension_vectors, sst,
  vgit_walls, wall_system, vgit_chambers, vgit_fan, git_equivalent, all_stability_parameters

import Pkg

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])

function _print_banner()
  printstyled(raw"""   ___"""; color=:red)
  printstyled(raw"""       _             """)
  println("  |")
  printstyled(raw"""  / _ \ """; color=:red)
  printstyled(raw"""_  _(_)_ _____ _ _ """)
  println("  |  Software package for quivers")
  printstyled(raw""" | (_) | """; color=:red)
  printstyled(raw"""|| | \ V / -_) '_|""")
  println("  |  and moduli of their representations")
  printstyled(raw"""  \__\_\\"""; color=:red)
  printstyled(raw"""\_,_|_|\_/\___|_|  """)
  println("  |")
  printstyled(raw"""       _____         _    """; color=:yellow)
  println("   |")
  printstyled(raw"""      |_   _|__  ___| |___"""; color=:yellow)
  println("   |  Manual: https://julia.quiver.tools")
  printstyled(raw"""        | |/ _ \/ _ \ (_-<"""; color=:yellow)
  println("   |  Version $(VERSION_NUMBER)")
  printstyled(raw"""        |_|\___/\___/_/__/"""; color=:yellow)
  return println("   |")
end

function __init__()
  if displaysize(stdout)[2] >= 80
    _print_banner()
  end

  return nothing
end

#######################################################
# Include all the submodules
#######################################################

include("Types.jl")
include("Quivers.jl")
include("Misc.jl")
include("Stability.jl")
include("RepresentationTheory.jl")
include("Constructors.jl")
include("Moduli.jl")
include("Hodge.jl")
include("Chow.jl")
include("Teleman.jl")
include("Bundles.jl")
include("WallsAndChambers.jl")

# Warm the JIT for the shared Chow/Hodge computation path so the user's first
# invariant computation is near-instant. Compilation is input-independent, so a
# single small example caches nearly all of it (see benchmark/ notes).
PrecompileTools.@compile_workload begin
  Q = kronecker_quiver(3)
  M = QuiverModuliSpace(Q, [2, 3])
  hodge_diamond(M)
  chow_ring(M)
  chern_numbers(M; unsafe=true)
end

######################
# end of QuiverTools
######################
end
