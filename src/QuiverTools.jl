module QuiverTools

using Pkg
using StaticArrays

using Memoization: Memoization
using IterTools: IterTools
using LinearAlgebraX: LinearAlgebraX
using Singular: Singular
using AbstractAlgebra: AbstractAlgebra
using Nemo: Nemo

import Base: show, ==, hash, getindex, length, iterate, keys, *, +, -, ^
import Memoization: @memoize, empty_all_caches!, empty_cache!
import IterTools: subsets
import LinearAlgebraX: rankx
import Singular: polynomial_ring, degree, coeff, constant_coefficient, AlgebraHomomorphism,
  preimage, Ideal, quotient_ideal, QuotientRing, fraction_field, std, gens, base_ring
import Combinatorics: combinations, with_replacement_combinations, partitions

# optional dependencies

try
  using SchubertPolynomials: SchubertPolynomials
  import SchubertPolynomials: xy_ring, schub_poly
catch
  @warn "SchubertPolynomials not found. Chow ring default functionnality may be slower.
  Solve this by `using Pkg; Pkg.add(url=\"https://github.com/pseudoeffective/SchubertPolynomials.jl\")`,
  then recompile QuiverTools.jl with `using Pkg; Pkg.build(\"QuiverTools\")."
end

# Types
export Quiver, HNType, LunaType, QuiverModuli, QuiverModuliSpace, QuiverModuliStack, Bundle

# Quivers
export n_vertices,
  n_arrows, arrows, indegree, outdegree, is_acyclic, is_connected, is_sink, is_source,
  underlying_graph, first_hochschild_cohomology

# Constructors
export kronecker_quiver, loop_quiver, subspace_quiver, three_vertex_quiver, cyclic_quiver,
  bipartite_quiver, opposite_quiver, double_quiver, dynkin_quiver
export kronecker_moduli, subspace_quiver_moduli

# Stability
export canonical_stability, is_coprime, slope
export all_hn_types,
  is_hn_type, has_semistables, has_stables, codimension_hn_stratum, is_amply_stable
export is_generic_subdimension_vector, all_generic_subdimension_vectors

# Representation theory
export euler_form, euler_matrix, is_schur_root, is_real_root, is_imaginary_root,
  is_isotropic_root,
  generic_ext, generic_hom, canonical_decomposition, in_fundamental_domain

# Moduli
export all_luna_types, is_luna_type, dimension_of_luna_stratum
export is_nonempty, codimension_unstable_locus, dimension, is_smooth,
  is_projective, semistable_equals_stable, semisimple_moduli_space

# Hodge
export hodge_diamond, hodge_polynomial, picard_rank, index, betti_numbers

# Chow
export chow_ring, motive, index, betti_numbers, poincare_polynomial, is_smooth,
  is_projective,
  semisimple_moduli_space, point_class, todd_class, chern_class_line_bundle,
  chern_character_line_bundle, total_chern_class_universal, integral

# Teleman
export teleman_bounds, weights_hn_type, weights_universal_bundle, weights_canonical_bundle,
  all_weights_endomorphisms_universal_bundle, weights_endomorphisms_universal_bundles,
  does_rigidity_inequality_hold, set_teleman_weights!

# Bundles
export chern_character, chern_class, chern_classes, dual, exterior_power, symmetric_power,
  det, line_bundle, canonical_bundle, universal_bundle, degree, rank, teleman_weights,
  structure_sheaf

# TODO add missing doctests across codebase
# TODO add safety checks everywhere in the codebase

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
include("Stability.jl")
include("RepresentationTheory.jl")
include("Misc.jl")
include("Constructors.jl")
include("Moduli.jl")
include("Hodge.jl")
include("Chow.jl")
include("Teleman.jl")
include("Bundles.jl")

######################
# end of QuiverTools
######################
end
