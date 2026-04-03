##############################################
# Quiver-specific helpers for AbstractBundle.
##############################################

"""
    chern_classes(F::Oscar.AbstractBundle)

Return the Chern classes of `F` as a dictionary keyed by the codimension.
"""
function chern_classes(F::Oscar.AbstractBundle)
  X = parent(F)
  return Dict(i => Oscar.chern_class(F, i) for i in 0:Oscar.dim(X))
end

"""
    teleman_weights(F::Oscar.AbstractBundle)

Return the Teleman weights stored on `F`.
"""
function teleman_weights(F::Oscar.AbstractBundle)
  weights = Oscar.get_attribute(F, :teleman_weights, nothing)
  isnothing(weights) && throw(ArgumentError("Bundle has no Teleman weights."))
  return weights
end

"""
    structure_sheaf(M::QuiverModuliSpace)

Return the structure sheaf of `M` as an `Oscar.AbstractBundle`.
"""
function structure_sheaf(M::QuiverModuliSpace)
  F = Oscar.trivial_line_bundle(chow_ring(M))
  return __attach_structure_sheaf_weights!(M, F)
end

function __line_bundle_class(
  M::QuiverModuliSpace,
  eta::AbstractVector{Int};
  unsafe::Bool=false,
)
  eta' * M.d != 0 && throw(ArgumentError("$(collect(eta)) is not a linearization."))

  X = chow_ring(M; unsafe=unsafe)
  A = Oscar.chow_ring(X)
  bundles = Oscar.tautological_bundles(X)
  return Oscar.simplify(
    -sum((eta[i] * Oscar.chern_class(bundles[i], 1) for i in support(M.d)); init=A(0))
  )
end

"""
    line_bundle(M::QuiverModuliSpace, eta::AbstractVector{Int}; unsafe::Bool=false, teleman::Bool=true)

Return the quiver line bundle defined by the linearization `eta` as an
`Oscar.AbstractBundle`.
"""
function line_bundle(
  M::QuiverModuliSpace,
  eta::AbstractVector{Int};
  unsafe::Bool=false,
  teleman::Bool=true,
)
  F = Oscar.line_bundle(
    chow_ring(M; unsafe=unsafe), __line_bundle_class(M, eta; unsafe=unsafe)
  )
  teleman && __attach_line_bundle_weights!(M, eta, F)
  return F
end

"""
    canonical_bundle(M::QuiverModuliSpace; teleman::Bool=true, unsafe::Bool=false)

Return the canonical bundle of `M` as an `Oscar.AbstractBundle`.
"""
function canonical_bundle(
  M::QuiverModuliSpace;
  teleman::Bool=true,
  unsafe::Bool=false,
)
  __validate_chow_inputs(M.Q, M.d, M.theta; unsafe=unsafe)
  F = line_bundle(M, -canonical_stability(M.Q, M.d); unsafe=true, teleman=false)
  teleman && __attach_canonical_weights!(M, F)
  return F
end

"""
    universal_bundle(M::QuiverModuliSpace, i::Int; teleman::Bool=true, unsafe::Bool=false)

Return the `i`-th tautological bundle of `M` as an `Oscar.AbstractBundle`.
"""
function universal_bundle(
  M::QuiverModuliSpace,
  i::Int;
  teleman::Bool=true,
  unsafe::Bool=false,
)
  X = chow_ring(M; unsafe=unsafe)
  1 <= i <= length(Oscar.tautological_bundles(X)) ||
    throw(BoundsError(Oscar.tautological_bundles(X), i))

  F = Oscar.tautological_bundles(X)[i]
  return teleman ? F : __copy_bundle(F)
end

"""
    degree(F::Oscar.AbstractBundle; unsafe::Bool=false)

Return the degree of `F`, defined as the top self-intersection of its first
Chern class.
"""
function degree(F::Oscar.AbstractBundle; unsafe::Bool=false)
  X = parent(F)
  return Oscar.integral(Oscar.chern_class(F, 1)^Oscar.dim(X))
end
