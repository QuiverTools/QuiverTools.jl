using Test, QuiverTools, Documenter
# To ensure that the version of the tests being run
# is the latest committed one,
# the local installation of QuiverTools should be the one developed at path="../."

# Tests that have mathematical significance
# should be in the documentation doctests.

@info "Almost all the tests are in the documentation."

# `import Oscar` (not `using`) loads Oscar so the walls-and-chambers / VGIT extension
# activates and its doctests resolve `Oscar.*`, without pulling Oscar's exports into
# scope (which would clash with QuiverTools names such as `index`, `todd_class`, ...).
DocMeta.setdocmeta!(QuiverTools, :DocTestSetup, :(using QuiverTools; import Oscar))
doctest(QuiverTools; manual=false, testset="Doctests")

@testset "strict sst" begin
  # proper-semistability
  Q = kronecker_quiver(2)
  d = [2, 2]
  theta = [1, -1]

  @test has_semistables(Q, d, theta) == true
  @test has_stables(Q, d, theta) == false
  Q = kronecker_quiver(3)
  @test has_stables(Q, [1, 0], [0, -1]) == true
  @test has_stables(Q, [0, 1], [1, 0]) == true
  @test has_stables(Q, [3, 0], [0, -1]) == false
  @test has_stables(Q, [0, 3], [1, 0]) == false
  @test has_semistables(Q, [3, 0], [0, -1]) == true
  @test has_semistables(Q, [0, 3], [1, 0]) == true
end;

@testset "HN types" begin
  # all_HN_types()

  Q = three_vertex_quiver(3, 4, 5)
  d = [3, 5, 7]
  theta = [43, 26, -37]

  # 3vertexquiver-3-5-7-canonical.txt
  expected = "" #has to be initialised outside of the open file
  open(joinpath(@__DIR__, "3vertexquiver-3-5-7-canonical.txt"), "r") do file
    expected = readline(file)
  end

  @test string(all_hn_types(Q, d, theta; ordered=true)) == expected
end;

@testset "Constructors" begin
  # equivalences from the Sage docstrings; == compares adjacency matrices only
  @test thickened_subspace_quiver(2, 6) == three_vertex_quiver(0, 6, 6)
  @test generalized_subspace_quiver(2, [2, 3]) == three_vertex_quiver(0, 2, 3)
  @test generalized_subspace_quiver(3, [1, 1, 1]) == subspace_quiver(3)
  @test jordan_quiver() == loop_quiver(1)

  U = disjoint_union(kronecker_quiver(3), kronecker_quiver(4))
  @test n_vertices(U) == 4
  @test n_arrows(U) == 7
end;

@testset "string constructor" begin
  # chain semantics: a run of r hyphens is r arrows, chains are read left to right.
  # integer tokens index vertices directly (so "1-6" is an arrow 1 -> 6); non-integer
  # tokens are labels numbered 1, ..., n in order of first appearance.
  @test Quiver("a---b") == kronecker_quiver(3)
  @test Quiver("1--2-3") == Quiver([0 2 0; 0 0 1; 0 0 0])
  @test Quiver("a--b-3,a---3,3-a") == Quiver([0 2 3; 0 0 1; 1 0 0])
  @test Quiver("1--2,1---3,2----3") == Quiver([0 2 3; 0 0 4; 0 0 0])

  # integer tokens fix the vertex indices; when they are exactly 1..n (in any order)
  # they are used as-is, with no warning (issue #28)
  @test (@test_nowarn Quiver("2-1")) == Quiver([0 0; 1 0])
  @test (@test_nowarn Quiver("1-6,2-6,3-6,4-6,5-6,7--1")) == Quiver(
    [
      0 0 0 0 0 1 0
      0 0 0 0 0 1 0
      0 0 0 0 0 1 0
      0 0 0 0 0 1 0
      0 0 0 0 0 1 0
      0 0 0 0 0 0 0
      2 0 0 0 0 0 0
    ],
  )
  # integer labels that are not 1..n are relabeled to 1..n (kth smallest -> vertex k),
  # with a warning: "2--3" and "1--3" both become "1--2"
  @test (@test_logs (:warn,) Quiver("2--3")) == Quiver([0 2; 0 0])
  @test (@test_logs (:warn,) Quiver("1--3")) == Quiver([0 2; 0 0])
  @test_throws ArgumentError Quiver("0-1")
end;

@testset "Luna types and local quivers" begin
  # is_luna_type must weight each key by the sum of its multiplicities (regression:
  # it previously summed the keys unweighted, wrongly rejecting any type with a
  # multiplicity list summing to more than one).
  Q = kronecker_quiver(3)
  M = QuiverModuliSpace(Q, [3, 3])
  for tau in all_luna_types(M)
    @test is_luna_type(M, tau.data)
  end
  @test is_luna_type(M, Dict([1, 1] => [3]))               # 3*[1,1] = [3,3]
  @test is_luna_type(M, Dict([1, 1] => [1], [2, 2] => [1])) # [1,1] + [2,2] = [3,3]
  @test !is_luna_type(M, Dict([1, 1] => [2]))              # 2*[1,1] = [2,2] != [3,3]

  # The encoding requires nonzero dimension vectors, nonempty lists of positive
  # multiplicities, and stable (not merely semistable) summands.
  @test !is_luna_type(M, Dict([1, 1] => Int[]))
  @test !is_luna_type(M, Dict([1, 1] => [-1], [2, 2] => [2]))
  @test !is_luna_type(M, Dict([1, 1, 0] => [3]))
  @test !is_luna_type(QuiverModuliSpace(kronecker_quiver(2), [2, 2]), Dict([2, 2] => [1]))

  # A rigid stable summand can occur with higher multiplicity, but there cannot be
  # two distinct stable summands of that dimension vector.
  R = QuiverModuliSpace(Q, [2, 0])
  @test is_luna_type(R, Dict([1, 0] => [2]))
  @test !is_luna_type(R, Dict([1, 0] => [1, 1]))
  @test_throws DomainError dimension_of_luna_stratum(R, Dict([1, 0] => [1, 1]))

  # Local quiver at a stable point is the g-loop quiver on one vertex with
  # g = 1 - <d,d> = dim M^s. For the 3-Kronecker quiver and d = (2,2) this is g = 5.
  # (This is the value from the definition in MR1972892; it intentionally differs from
  # QuiverTools/Sage, which returns 4 via general_ext and undercounts the diagonal.)
  X = QuiverModuliSpace(Q, [2, 2])
  loc = QuiverTools.local_quiver_setting(X, Dict([2, 2] => [1]))
  @test propertynames(loc) == (:Q, :d, :summands)
  @test loc.d == [1]
  @test Matrix(loc.Q.adjacency) == fill(5, 1, 1)

  # Luna strata of M(2d) ≅ P^2 for the subspace quiver Q^(4) = affine D4 with
  # d = (1,1,1,1;2). There are five polystable types; we check their local quivers
  # and local dimension vectors against a hand computation.
  S = subspace_quiver(4)
  d = [1, 1, 1, 1, 2]
  N = QuiverModuliSpace(S, 2 .* d)
  # fingerprint invariant under relabeling the local quiver's vertices: the sorted
  # local dimension vector, the sorted loop counts, the sorted adjacency entries, and
  # whether the local quiver is symmetric (i.e. only loops and 2-cycles).
  function fingerprint(tau)
    s = QuiverTools.local_quiver_setting(N, tau)
    A = Matrix(s.Q.adjacency)
    loops = [A[i, i] for i in 1:size(A, 1)]
    (sort(s.d), sort(loops), sort(vec(A)), A == permutedims(A))
  end
  eK, eKb, eL, eLb = [1, 1, 0, 0, 1], [0, 0, 1, 1, 1], [1, 0, 1, 0, 1], [0, 1, 0, 1, 1]
  # ξ1 = (d, d): two vertices, a loop on each, no arrows between them
  @test fingerprint(Dict(d => [1, 1])) == ([1, 1], [1, 1], [0, 0, 1, 1], true)
  # ξ2 = (d, e_K, e_Kbar): a loop on the d-vertex and a 2-cycle on (e_K, e_Kbar)
  @test fingerprint(Dict(d => [1], eK => [1], eKb => [1])) ==
    ([1, 1, 1], [0, 0, 1], sort([1, 1, 1, 0, 0, 0, 0, 0, 0]), true)
  # ξ3 = (e_K, e_Kbar, e_L, e_Lbar) with |K∩L| = 1: two disjoint 2-cycles
  @test fingerprint(Dict(eK => [1], eKb => [1], eL => [1], eLb => [1])) ==
    ([1, 1, 1, 1], [0, 0, 0, 0], sort(vcat(fill(1, 4), fill(0, 12))), true)
  # ξ4 = (d^2): one vertex with a single loop and local dimension 2
  @test fingerprint(Dict(d => [2])) == ([2], [1], [1], true)
  # ξ5 = (e_K^2, e_Kbar^2): two vertices, a 2-cycle, local dimension (2, 2)
  @test fingerprint(Dict(eK => [2], eKb => [2])) == ([2, 2], [0, 0], [0, 0, 1, 1], true)
end;

@testset "Betti numbers" begin
  # betti_numbers is the coefficient vector of the Poincaré polynomial, indexed
  # by ascending cohomological degree, always a Vector{Int} (regression: for
  # non-palindromic Poincaré polynomials, which occur for quivers with oriented
  # cycles, the coefficients were listed by descending degree and the padding to
  # length 2 dim + 1 used Float64 zeros; a monomial Poincaré polynomial crashed).
  # The round-trip quiver with d = (1, 1) has moduli space A^1, so P = L.
  Q = Quiver([0 1; 1 0])
  M = QuiverModuliSpace(Q, [1, 1], [1, -1])
  @test betti_numbers(M) isa Vector{Int}
  @test betti_numbers(M) == [0, 0, 1]

  # the smooth projective (palindromic) case is unchanged: our favourite 6-fold
  M = QuiverModuliSpace(kronecker_quiver(3), [2, 3])
  @test betti_numbers(M) == [1, 0, 1, 0, 3, 0, 3, 0, 3, 0, 1, 0, 1]
end;

@testset "motive/poincare sign" begin
  # The two branches of motive must agree in sign: the HN recursion returns the
  # honest stack motive [M]/(L-1), and the trivial-stability branch returns the
  # full-stack motive; poincare_polynomial recovers [M] = (L-1)*[stack motive].
  # (regression #36: for a quiver with loops the trivial branch was consumed with
  # the opposite sign, giving e.g. P = -L for the affine line.)
  # poincare_polynomial returns an spoly and motive an n_transExt, so compare the
  # rendered strings, as the doctests do.

  # loop_quiver(g), d = [1]: moduli space is A^g, so P = L^g and motive = L^g/(L-1)
  for (g, p) in ((1, "L"), (2, "L^2"), (3, "L^3"))
    M = QuiverModuliSpace(loop_quiver(g), [1], [0])
    @test string(poincare_polynomial(M)) == p
    @test string(motive(QuiverModuliStack(loop_quiver(g), [1], [0], "stable"))) ==
      "$p//(L - 1)"
  end

  # recursion path is unchanged: the space motive of the 6-fold stays positive
  M = QuiverModuliSpace(kronecker_quiver(3), [2, 3])
  @test string(poincare_polynomial(M)) ==
    "L^6 + L^5 + 3*L^4 + 3*L^3 + 3*L^2 + L + 1"
end;

@testset "Bocklandt reduction" begin
  # The public quiver-setting routines enforce the dimension-vector contract.
  for f in (bocklandt_reduction, is_coregular, is_cofree)
    @test_throws ArgumentError f(jordan_quiver(1), [1, 1])
    @test_throws ArgumentError f(jordan_quiver(1), [-1])
  end

  # The dimension API rejects disconnected quivers with the documented exception.
  disconnected = disjoint_union(kronecker_quiver(1), kronecker_quiver(1))
  @test_throws ArgumentError dimension(QuiverModuliSpace(disconnected, [1, 1, 1, 1]))

  # invariants of pairs of 2x2 matrices form a polynomial ring, of 3x3 they do not,
  # and neither do those of triples of 2x2 matrices; a single matrix always does
  @test is_coregular(jordan_quiver(2), [2])
  @test !is_coregular(jordan_quiver(2), [3])
  @test !is_coregular(jordan_quiver(3), [2])
  @test all(is_coregular(jordan_quiver(1), [n]) for n in 1:5)

  # for acyclic quivers the quotient variety is a point
  @test is_coregular(kronecker_quiver(3), [2, 3])
  @test is_coregular(subspace_quiver(4), [1, 1, 1, 1, 2])

  # settings I, II and IV of [Theorem 4.4, MR1929191] are coregular
  @test is_coregular(Quiver("1-2, 2-1"), [4, 5])              # I
  @test is_coregular(Quiver("1--2, 2--1"), [1, 2])            # II with k = 2 <= n = 2
  @test !is_coregular(Quiver("1--2, 2--1"), [1, 1])           # II fails for k = 2 > n = 1
  @test is_coregular(Quiver("1-2, 2-1, 2-3, 3-2"), [3, 2, 3]) # IV

  # the reduction combines R_III, R_I and R_II to a lone vertex of dimension 1
  setting = bocklandt_reduction(Quiver("1-2, 2-2, 2-1"), [1, 2])
  @test propertynames(setting) == (:Q, :d)
  @test n_vertices(setting.Q) == 1
  @test n_arrows(setting.Q) == 0
  @test setting.d == [1]

  # a reduced setting is returned unchanged
  setting = bocklandt_reduction(Quiver("1--2, 2--1"), [1, 1])
  @test Matrix(setting.Q.adjacency) == [0 2; 2 0]
  @test setting.d == [1, 1]

  # vertices of dimension 0 and arrows between strongly connected components are dropped
  @test bocklandt_reduction(kronecker_quiver(3), [2, 0]).d == [2]
  @test is_coregular(Quiver("1-1, 1-2, 2-2"), [2, 2])
  @test is_coregular(kronecker_quiver(3), [0, 0])

  # smoothness of moduli spaces with properly semistable representations: for the
  # 2-Kronecker quiver and d = (2, 2) one gets P^2, and for the 3-Kronecker quiver
  # both d = (2, 2) and d = (2, 4) give P^5: the deepest local quiver setting is two
  # loops on a vertex of dimension 2, the reduced coregular setting C1 of [MR1929191];
  # for d = (3, 3) that setting has dimension 3 instead, so the space is singular
  @test is_smooth(QuiverModuliSpace(kronecker_quiver(2), [2, 2]))
  @test is_smooth(QuiverModuliSpace(kronecker_quiver(3), [2, 2]))
  @test is_smooth(QuiverModuliSpace(kronecker_quiver(3), [2, 4]))
  @test !is_smooth(QuiverModuliSpace(kronecker_quiver(3), [3, 3]))

  # the 6-subspace quiver with d = (1^5, 2; 3) and stability parameters on a wall:
  # for theta = (1^5, 2; -3) the moduli space is accidentally isomorphic to Gr(2, 4),
  # hence smooth despite the eleven Luna strata, whereas for theta = (2^5, 1; -4)
  # there are ten isolated singular points, one for each two-element subset of the
  # five thin subspace vertices
  S = subspace_quiver(6)
  d = [1, 1, 1, 1, 1, 2, 3]
  @test is_smooth(QuiverModuliSpace(S, d, [1, 1, 1, 1, 1, 2, -3]))
  @test !is_smooth(QuiverModuliSpace(S, d, [2, 2, 2, 2, 2, 1, -4]))
  # for d = (1^4, 2^2; 3) the analogous first wall crossing has smooth target too
  @test is_smooth(QuiverModuliSpace(S, [1, 1, 1, 1, 2, 2, 3], [2, 2, 2, 2, 1, 1, -4]))

  # the Segre cubic threefold, as the moduli space for the 6-subspace quiver with
  # d = (1^6; 2) and canonical stability: it has ten nodes, one for each splitting
  # of the six thin subspace vertices into complementary triples
  @test !is_smooth(QuiverModuliSpace(S, [1, 1, 1, 1, 1, 1, 2]))

  # codimension of the singular locus: the ten nodes of the Segre cubic; for the
  # 3-Kronecker quiver and d = (2, 2) the properly semistable locus is non-empty
  # while the singular locus is empty, and for d = (3, 3) the largest singular
  # Luna stratum has codimension 3 in the 10-dimensional moduli space
  @test codimension_singular_locus(QuiverModuliSpace(S, [1, 1, 1, 1, 1, 1, 2])) == 3
  @test codimension_singular_locus(QuiverModuliSpace(kronecker_quiver(3), [2, 2])) == Inf
  @test codimension_singular_locus(QuiverModuliSpace(kronecker_quiver(3), [3, 3])) == 3
end;

@testset "cofree quiver settings" begin
  # cyclic quiver settings and matrix invariants: pairs of 2x2 matrices are cofree,
  # pairs of 3x3 matrices and triples of 2x2 matrices are not; any number of loops on
  # a vertex of dimension 1 is cofree
  @test all(is_cofree(cyclic_quiver(n), fill(k, n)) for n in 1:3, k in 1:3)
  @test is_cofree(cyclic_quiver(3), [1, 2, 3])
  @test is_cofree(jordan_quiver(2), [2])
  @test !is_cofree(jordan_quiver(2), [3])
  @test !is_cofree(jordan_quiver(3), [2])
  @test is_cofree(jordan_quiver(3), [1])

  # acyclic settings are trivially cofree, as the invariants are constants
  @test is_cofree(kronecker_quiver(3), [2, 3])
  @test is_cofree(subspace_quiver(4), [1, 1, 1, 1, 2])

  # settings with all cycles through a vertex of dimension 1: the k arrows back and
  # forth give 2k - 1 as the bound on the other dimension
  @test is_cofree(Quiver("1-2, 2-1"), [1, 5])
  @test is_cofree(Quiver("1--2, 2--1"), [1, 3])
  @test !is_cofree(Quiver("1--2, 2--1"), [1, 2])

  # two cycles sharing a path: cofree iff exactly one shared dimension is 2 and the
  # others are at least 4, so coregularity does not suffice
  @test is_coregular(Quiver("1--2, 2-1"), [2, 2])
  @test !is_cofree(Quiver("1--2, 2-1"), [2, 2])
  @test !is_cofree(Quiver("1--2, 2-1"), [2, 3])
  @test is_cofree(Quiver("1--2, 2-1"), [2, 4])

  # two cycles sharing a path through a vertex of dimension 1: cofree iff the minimal
  # dimension along the big cycle is attained exactly once in the shared path, or not
  # there but exactly once in the other branch
  theta_quiver = Quiver("1-2, 2-3, 3-1, 2-4, 4-1")
  @test is_cofree(theta_quiver, [2, 3, 4, 1])
  @test is_cofree(theta_quiver, [3, 3, 2, 1])
  @test !is_cofree(theta_quiver, [2, 2, 3, 1])

  # wedging removes the vertex of dimension 3 on the path to the central vertex,
  # reducing to the setting [2, 3, 4, 1] above; with dimension 1 instead there are two
  # vertices of dimension 1 on a common cycle, which is never cofree
  wedged = Quiver("1-2, 2-3, 3-1, 2-5, 5-4, 4-1")
  @test is_cofree(wedged, [2, 3, 4, 1, 3])
  @test !is_cofree(wedged, [2, 3, 4, 1, 1])
end;

@testset "framed quiver moduli" begin
  # Quiver-level framing (Sage parity): the framing vertex is prepended (arrows i0 -> i),
  # the coframing vertex appended (arrows i -> i0). See arXiv:2607.12895.
  kro = kronecker_quiver(2)
  @test Matrix(framed_quiver(kro, [0, 1]).adjacency) == [0 0 1; 0 0 2; 0 0 0]
  @test Matrix(coframed_quiver(kro, [0, 1]).adjacency) == [0 2 0; 0 0 1; 0 0 0]
  @test_throws ArgumentError framed_quiver(kro, [1])

  # A framed quiver moduli space: `base` is the codomain of the projection p, `total_space`
  # the framed moduli space as an ordinary QuiverModuliSpace.
  X = FramedQuiverModuliSpace(kro, [2, 2]; n=[0, 2])
  @test base(X).d == [2, 2]
  @test framing_vector(X) == [0, 2]
  @test n_vertices(total_space(X).Q) == 3

  # local_quiver_setting is public API (issue #41): the local quiver and dimension vector of
  # a Luna stratum. For the 3-Kronecker quiver and the type [1,1] with multiplicity 3, the
  # local quiver is the two-loop quiver on one vertex.
  M = QuiverModuliSpace(kronecker_quiver(3), [3, 3])
  @test local_quiver_setting(M, Dict([1, 1] => [3])).d == [3]

  # Proposition 5: the etale-local model of the singularity along a Luna stratum is the local
  # quiver at the trivial stability parameter.
  L = local_structure(M, Dict([1, 1] => [3]))
  @test L isa QuiverModuliSpace
  @test L.d == [3]
  @test Matrix(L.Q.adjacency) == fill(2, 1, 1)   # one vertex, two loops
  @test L.theta == [0]

  # Proposition 6: fibres of the framed projection for the framed affine-D4 quiver, base
  # dimension 2*(1,1,1,1;2). This is the four-dimensional QM_1 of arXiv:2607.12895, whose
  # five fibre types are Lemmas 18-22. The fibre over each stratum is a NilpotentLocus in a
  # framed moduli of the local quiver, with local framing datum n_tau = sum_k (n.d_k) i_k and
  # framed stability [sum d_tau, -1, ..., -1] (framing vertex first). Local vertex order is
  # not canonical, so the framing datum is compared sorted. Rows: (Lemma, tau, d_tau, n_tau).
  delta = [1, 1, 1, 1, 2]
  XD4 = FramedQuiverModuliSpace(subspace_quiver(4), 2 .* delta; n=[0, 0, 0, 0, 1])
  @test dimension(total_space(XD4)) == 4
  eK, eKb = [1, 1, 0, 0, 1], [0, 0, 1, 1, 1]
  eL, eLb = [1, 0, 1, 0, 1], [0, 1, 0, 1, 1]
  types = [
    (18, Dict(delta => [1, 1]), [1, 1], [2, 2]),
    (19, Dict(delta => [1], eK => [1], eKb => [1]), [1, 1, 1], [1, 1, 2]),
    (20, Dict(eK => [1], eKb => [1], eL => [1], eLb => [1]), [1, 1, 1, 1], [1, 1, 1, 1]),
    (21, Dict(delta => [2]), [2], [2]),
    (22, Dict(eK => [2], eKb => [2]), [2, 2], [1, 1]),
  ]
  for (lemma, tau, dloc, nloc) in types
    F = fibre(XD4, tau)
    @test F isa NilpotentLocus
    A = ambient(F)
    @test A isa FramedQuiverModuliSpace
    @test sort(A.d) == dloc
    @test sort(framing_vector(A)) == nloc
    @test total_space(A).theta == [sum(dloc); fill(-1, length(dloc))]
  end
end;
