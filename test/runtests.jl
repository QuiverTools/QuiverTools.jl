using Test, QuiverTools, Documenter, Logging
# To ensure that the version of the tests being run
# is the latest committed one,
# the local installation of QuiverTools should be the one developed at path="../."

# Tests that have mathematical significance
# should be in the documentation doctests.

@testset "Oscar-independent VGIT predicates" begin
  # VGIT boolean queries must not load the optional Oscar extension.
  @test Base.get_extension(QuiverTools, :QuiverToolsOscarExt) === nothing

  Q = three_vertex_quiver(2, 3, 4)
  d = [1, 2, 2]
  theta1 = [2, -1 // 2, -1 // 2]
  theta2 = [2, 1 // 2, -3 // 2]
  @test !git_equivalent(Q, d, theta1, theta2)
  @test git_equivalent(Q, d, theta1, 2 .* theta1)

  S = subspace_quiver(6)
  d = [1, 1, 1, 1, 1, 2, 3]
  theta = [3, 3, 3, 3, 3, 3, -7]
  thetabar = [1, 1, 1, 1, 1, 2, -3]
  @test QuiverTools.__is_vgit_chamber_parameter(S, d, theta)
  @test QuiverTools.__in_closure_of_vgit_chamber(S, d, theta, thetabar)
  @test !QuiverTools.__is_vgit_chamber_parameter(S, d, thetabar)

  @test Base.get_extension(QuiverTools, :QuiverToolsOscarExt) === nothing
end;
QuiverTools.empty_all_caches!()

# Oscar has to be imported here to not throw off the previous test set
import Oscar

@testset "VGIT predicates agree with Oscar " begin
  Q = three_vertex_quiver(2, 3, 4)
  d = [1, 2, 2]
  parameters = [
    [2, -1 // 2, -1 // 2],
    [2, 1 // 2, -3 // 2],
    [1, 3 // 2, -2],
    [0, 1, -1],
    [1, 0, -1 // 2],
    [0, 0, 0],
  ]

  cone = sst(Q, d)
  for theta in parameters
    @test QuiverTools.__in_semistable_cone(Q, d, theta) == (theta in cone)
  end

  for e in QuiverTools.all_subdimension_vectors(d; nonzero=true, strict=true)
    wall = reduce(Oscar.intersect, (sst(Q, e), sst(Q, d - e), cone))
    for theta in parameters
      @test QuiverTools.__in_vgit_wall(Q, d, e, theta) == (theta in wall)
    end
  end

  walls = vgit_walls(Q, d; top_dimension=false)
  function oscar_git_equivalent(theta1, theta2)
    theta1 == theta2 && return true
    line = Oscar.convex_hull([theta1, theta2])
    return all(walls) do wall
      Oscar.issubset(line, wall) || !Oscar.is_feasible(Oscar.intersect(line, wall))
    end
  end
  for theta1 in parameters, theta2 in parameters
    @test git_equivalent(Q, d, theta1, theta2) ==
      oscar_git_equivalent(theta1, theta2)
  end
end;
QuiverTools.empty_all_caches!()

@testset "vgit_chambers edge cases         " begin
  # a single vertex with d = [1] is a primitive Schur root and must not be rejected;
  # fixed in 15b4405
  @test length(vgit_chambers(Quiver(zeros(Int, 1, 1)), [1])) == 1

  # splitting along earlier walls cuts some of the smaller walls into several pieces,
  # each of which must still split the chamber it meets (this used to give 5 chambers);
  # fixed in 783271c
  @test length(vgit_chambers(Quiver("1-2,1-3,1--4,2--3,2-4,3-4"), [1, 1, 1, 1])) == 7

  # here one chamber is split into several pieces by fake walls; after merging two of
  # them, the fake-wall loop kept iterating over the old list of pieces, so a piece could
  # be merged a second time, giving overlapping chambers (11 instead of 9);
  # fixed in 24b4a9b
  Q = Quiver("1--4,2-3,3-4,5--3")
  @test length(vgit_chambers(Q, [1, 1, 2, 2, 2])) == 9
end;
QuiverTools.empty_all_caches!()

@testset "vgit_fan with lineality spaces   " begin
  # vgit_fan used to reject chambers with a nontrivial lineality space

  # d has a zero entry, so the chamber contains a line
  F = vgit_fan(Quiver([0 3 0; 0 0 0; 0 0 0]), [2, 3, 0])
  @test (Oscar.dim(F), Oscar.lineality_dim(F), Oscar.n_maximal_cones(F)) == (2, 1, 1)

  # there are no walls and the only chamber, sst(Q, d), is a linear space
  F = vgit_fan(subspace_quiver(3), [0, 0, 0, 1])
  @test (Oscar.dim(F), Oscar.lineality_dim(F), Oscar.n_maximal_cones(F)) == (3, 3, 1)
end;
QuiverTools.empty_all_caches!()

@testset "strict sst                       " begin
  # proper-semistability
  Q = kronecker_quiver(2)
  d = [2, 2]
  theta = [1, -1]

  @test has_semistable_representation(Q, d, theta) == true
  @test has_stable_representation(Q, d, theta) == false
  Q = kronecker_quiver(3)
  @test has_stable_representation(Q, [1, 0], [0, -1]) == true
  @test has_stable_representation(Q, [0, 1], [1, 0]) == true
  @test has_stable_representation(Q, [3, 0], [0, -1]) == false
  @test has_stable_representation(Q, [0, 3], [1, 0]) == false
  @test has_semistable_representation(Q, [3, 0], [0, -1]) == true
  @test has_semistable_representation(Q, [0, 3], [1, 0]) == true
end;
QuiverTools.empty_all_caches!()

@testset "HN types                         " begin
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
QuiverTools.empty_all_caches!()

@testset "Constructors                     " begin
  # equivalences from the Sage docstrings; == compares adjacency matrices only
  @test thickened_subspace_quiver(2, 6) == three_vertex_quiver(0, 6, 6)
  @test generalized_subspace_quiver(2, [2, 3]) == three_vertex_quiver(0, 2, 3)
  @test generalized_subspace_quiver(3, [1, 1, 1]) == subspace_quiver(3)
  @test jordan_quiver() == loop_quiver(1)

  U = disjoint_union(kronecker_quiver(3), kronecker_quiver(4))
  @test n_vertices(U) == 4
  @test n_arrows(U) == 7
end;
QuiverTools.empty_all_caches!()

@testset "string constructor               " begin
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
QuiverTools.empty_all_caches!()

@testset "Luna types and local quivers     " begin
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
QuiverTools.empty_all_caches!()

@testset "Betti numbers                    " begin
  # betti_numbers is the coefficient vector of the Poincaré polynomial, indexed
  # by ascending cohomological degree, always a Vector{Int} (regression: for
  # non-palindromic Poincaré polynomials, which occur for quivers with oriented
  # cycles, the coefficients were listed by descending degree and the padding to
  # length 2 dim + 1 used Float64 zeros; a monomial Poincaré polynomial crashed).
  # The round-trip quiver with d = (1, 1) has moduli space A^1, so P = L.
  Q = Quiver([0 1; 1 0])
  M = QuiverModuliSpace(Q, [1, 1], [1, -1])
  # the quiver has a cycle, so this is cohomology with compact support; only with
  # `verbose = true` a warning says so
  betti = @test_logs betti_numbers(M)
  @test betti isa Vector{Int}
  @test betti == [0, 0, 1]
  @test_logs (:warn, r"not acyclic") betti_numbers(M; verbose=true)

  # the smooth projective (palindromic) case is unchanged: our favourite 6-fold
  M = QuiverModuliSpace(kronecker_quiver(3), [2, 3])
  @test betti_numbers(M) == [1, 0, 1, 0, 3, 0, 3, 0, 3, 0, 1, 0, 1]
end;
QuiverTools.empty_all_caches!()

@testset "motive/poincare sign             " begin
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
    @test string(@test_logs poincare_polynomial(M)) == p
    @test string(motive(QuiverModuliStack(loop_quiver(g), [1], [0], "stable"))) ==
      "$p//(L - 1)"
  end

  # recursion path is unchanged: the space motive of the 6-fold stays positive
  M = QuiverModuliSpace(kronecker_quiver(3), [2, 3])
  @test string(poincare_polynomial(M)) ==
    "L^6 + L^5 + 3*L^4 + 3*L^3 + 3*L^2 + L + 1"
end;
QuiverTools.empty_all_caches!()

@testset "theta outside d^perp             " begin
  # theta and theta + c * (1, ..., 1) define the same stability for denom = sum, so a
  # stability parameter need not lie in d^perp;

  K2, K3 = kronecker_quiver(2), kronecker_quiver(3)

  # [2, 1] = [1, -1] + 3/2 * (1, 1) for K2 and d = [2, 2]: properly semistable
  @test has_properly_semistable_representation(K2, [2, 2], [2, 1])
  @test !semistable_equals_stable(QuiverModuliSpace(K2, [2, 2], [2, 1]))
  @test_throws ArgumentError poincare_polynomial(K2, [2, 2], [2, 1])

  # [4, -1] = [3, -2] + (1, 1) for K3 and d = [2, 3]: our favourite 6-fold
  M, N = QuiverModuliSpace(K3, [2, 3], [4, -1]), QuiverModuliSpace(K3, [2, 3])
  @test !has_properly_semistable_representation(K3, [2, 3], [4, -1])
  @test semistable_equals_stable(M)
  @test poincare_polynomial(M) == poincare_polynomial(N)
  @test hodge_diamond(M) == hodge_diamond(N)
  @test betti_numbers(M) == betti_numbers(N)
  @test picard_rank(M) == 1
  @test index(M) == 3
end;
QuiverTools.empty_all_caches!()

@testset "intersection cohomology          " begin
  # intersection cohomology agrees with ordinary cohomology on smooth moduli spaces, so
  # route a few of those through the Meinhardt--Reineke branch on purpose
  for (Q, d, theta) in [
    (kronecker_quiver(2), [1, 1], [1, -1]),
    (kronecker_quiver(5), [1, 1], nothing),
    (kronecker_quiver(4), [1, 2], [2, -1]),
    (kronecker_quiver(5), [2, 3], nothing),
    (Quiver([0 1 1; 0 0 2; 0 0 0]), [1, 1, 1], nothing),
  ]
    M = if isnothing(theta)
      QuiverModuliSpace(Q, d)
    else
      QuiverModuliSpace(Q, d, theta)
    end
    @test intersection_betti_numbers(M) == betti_numbers(M)
  end

  # reflection functors identify M(a, b) with M(b, m*b - a) for the m-Kronecker quiver,
  # on dimension vectors whose moduli space is singular
  for (m, d, e) in [
    (3, [2, 2], [2, 4]),
    (3, [3, 3], [3, 6]),
    (3, [2, 3], [3, 7]),
    (4, [3, 3], [3, 9]),
    (5, [2, 2], [2, 8]),
  ]
    Q = kronecker_quiver(m)
    @test intersection_poincare_polynomial(QuiverModuliSpace(Q, d)) ==
      intersection_poincare_polynomial(QuiverModuliSpace(Q, e))
  end

  # intersection cohomology of a projective variety satisfies Poincaré duality
  for (Q, d) in [
    (kronecker_quiver(3), [2, 2]),
    (kronecker_quiver(4), [2, 4]),
    (Quiver([0 1 1; 0 0 2; 0 0 0]), [2, 2, 2]),
  ]
    M = QuiverModuliSpace(Q, d)
    betti = intersection_betti_numbers(M)
    @test betti == reverse(betti)
    @test all(b >= 0 for b in betti)
    @test length(betti) == 2 * dimension(M) + 1
  end

  # a quiver need not be acyclic, but then the moduli space is affine rather than
  # projective and the answer is compactly supported. These are the spaces of matrix
  # invariants, whose intersection cohomology is Theorem 8.2 of [MR4000572]; for `d = 1`
  # the moduli space is A^m, and for m = d = 2 it is A^5 because the five traces and
  # determinants are independent.
  for (m, d, expected) in [
    (2, 1, [2]),
    (3, 1, [3]),
    (2, 2, [5]),
    (2, 3, [10]),
    (2, 4, [17, 15]),
    (3, 2, [9]),
    (3, 3, [19, 17, 16]),
    (4, 2, [13, 11]),
  ]
    M = QuiverModuliSpace(loop_quiver(m), [d])
    @test dimension(M) == (m - 1) * d^2 + 1
    betti = intersection_betti_numbers(M)
    @test findall(!iszero, betti) .- 1 == sort(2 .* expected)
    @test all(isone, betti[findall(!iszero, betti)])
  end

  # the Donaldson--Thomas invariant vanishes when nothing of dimension vector `d` is
  # stable, which is an exact cancellation in the plethystic logarithm
  for n in 2:4
    M = QuiverModuliSpace(kronecker_quiver(2), [n, n])
    @test_throws ArgumentError intersection_poincare_polynomial(M)
  end

  # the theorem needs the stability parameter to be generic for the slope of `d`
  M = QuiverModuliSpace(Quiver([0 1 1 0; 0 0 1 0; 0 0 0 1; 0 0 0 0]), [3, 3, 4, 1])
  @test_throws ArgumentError intersection_poincare_polynomial(M)

  # and it describes the semistable moduli space
  M = QuiverModuliSpace(kronecker_quiver(3), [2, 3], [3, -2], "stable")
  @test_throws ArgumentError intersection_poincare_polynomial(M)
end;
QuiverTools.empty_all_caches!()

@testset "poincare via intersection        " begin
  # with properly semistable representations poincare_polynomial falls back on
  # intersection cohomology, which is ordinary cohomology when the moduli space is smooth
  K2, K3 = kronecker_quiver(2), kronecker_quiver(3)

  # K3 with d = [2, 2] is P^5
  M = QuiverModuliSpace(K3, [2, 2])
  identity = QuiverTools.identity_matrix(6)
  @test string(poincare_polynomial(M)) == "L^5 + L^4 + L^3 + L^2 + L + 1"
  @test betti_numbers(M) == [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
  @test picard_rank(M) == 1
  @test hodge_diamond(M) == identity
  @test hodge_diamond(K3, [2, 2]) == identity
  @test hodge_polynomial(K3, [2, 2]) == hodge_polynomial(M)

  # singular, and smooth (P^2) without stable representations: still refused
  @test_throws ArgumentError poincare_polynomial(QuiverModuliSpace(K3, [3, 3]))
  @test_throws ArgumentError poincare_polynomial(QuiverModuliSpace(K2, [2, 2]))

  # the stable moduli space is not the semistable one, so the fallback does not apply
  S = QuiverModuliSpace(K3, [2, 2], canonical_stability(K3, [2, 2]), "stable")
  @test_throws ArgumentError poincare_polynomial(S)
  @test_throws ArgumentError hodge_polynomial(S)
  @test_throws ArgumentError hodge_diamond(S)
  @test_throws ArgumentError betti_numbers(S)
  # unless there are no properly semistable representations
  @test poincare_polynomial(QuiverModuliSpace(K3, [2, 3], [3, -2], "stable")) ==
    poincare_polynomial(QuiverModuliSpace(K3, [2, 3]))

  # two loops and d = [2]: matrix invariants of pairs of 2x2 matrices, which is A^5, so
  # cohomology with compact support, pointed out in verbose mode
  C = QuiverModuliSpace(loop_quiver(2), [2], [0])
  @test is_smooth(C)
  @test string(@test_logs (:warn, r"not acyclic") poincare_polynomial(C; verbose=true)) ==
    "L^5"
  @test string(@test_logs poincare_polynomial(C)) == "L^5"
end;
QuiverTools.empty_all_caches!()

@testset "to_dot                           " begin
  # the exact string for the happy path is locked by the doctests; here we check
  # arrow bookkeeping and the edge cases they don't show.

  # one `->` line per arrow, loops included
  count_arrows(Q) = count(l -> occursin("->", l), split(to_dot(Q), '\n'))
  @test count_arrows(kronecker_quiver(3)) == n_arrows(kronecker_quiver(3)) == 3
  @test count_arrows(three_vertex_quiver(1, 2, 3)) ==
    n_arrows(three_vertex_quiver(1, 2, 3)) == 6
  @test count_arrows(loop_quiver(4)) == n_arrows(loop_quiver(4)) == 4

  # self-loops are spread with compass ports rather than stacked
  @test occursin("1:n -> 1:n", to_dot(loop_quiver(2)))
  @test occursin("1:s -> 1:s", to_dot(loop_quiver(2)))

  # a quiver with no arrows still declares its (isolated) vertices, no crash
  isolated = to_dot(Quiver([0 0; 0 0]))
  @test count(l -> occursin("->", l), split(isolated, '\n')) == 0
  @test occursin("\n  1;\n", isolated) && occursin("\n  2;\n", isolated)

  # graph names are quoted DOT strings, so syntax characters and physical line
  # breaks must be escaped before interpolation; UTF-8 and tabs remain literal
  named = Quiver(zeros(Int, 1, 1), "α → β, 日本語: quote \" slash \\ line\nbreak\ttab")
  @test occursin(
    "label=\"α → β, 日本語: quote \\\" slash \\\\ line\\nbreak\ttab\";",
    to_dot(named),
  )
end;
QuiverTools.empty_all_caches!()

import GraphViz

@testset "GraphViz extension               " begin
  Q = Quiver([0 1; 0 0], "α → β, 日本語: quote \" slash \\ line\nbreak")
  svg = sprint(show, MIME("image/svg+xml"), Q)
  @test occursin("<svg", svg)
  @test occursin("α → β, 日本語", svg)

  extension_module = Base.get_extension(QuiverTools, :QuiverToolsGraphVizExt)
  @test !isnothing(extension_module)
  mktempdir() do directory
    path = joinpath(directory, "quiver.svg")
    @test extension_module._write_svg(Q, path) == path
    @test isfile(path)
    @test occursin("<svg", read(path, String))
  end
end;
QuiverTools.empty_all_caches!()

@testset "Bocklandt reduction              " begin
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
QuiverTools.empty_all_caches!()

@testset "cofree quiver settings           " begin
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
QuiverTools.empty_all_caches!()

@testset "framed quiver moduli             " begin
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
end;
QuiverTools.empty_all_caches!()

@testset "nullcones and defect             " begin
  # The public nullcone routines enforce the dimension-vector contract.
  for f in (dimension_nullcone, nullcone_motive, defect)
    @test_throws ArgumentError f(jordan_quiver(1), [1, 1])
    @test_throws ArgumentError f(jordan_quiver(1), [-1])
  end

  # nilpotent matrices, and pairs of matrices with a common complete flag
  @test dimension_nullcone(jordan_quiver(1), [3]) == 6
  @test dimension_nullcone(jordan_quiver(2), [2]) == 3
  # acyclic settings have no invariants, so the nullcone is everything
  @test dimension_nullcone(kronecker_quiver(3), [2, 3]) == 18
  # the motivic recursion reaches larger settings, and for symmetric quivers it
  # matches the closed formula of [Remark 3.6, doi:10.3842/SIGMA.2026.020]
  @test dimension_nullcone(jordan_quiver(2), [4]) == 18
  @test dimension_nullcone(jordan_quiver(3), [2]) == 4
  @test dimension_nullcone(Quiver("1-1, 1-2, 2-1"), [2, 3]) == 11
  # the motive of pairs of nilpotent 3x3 matrices is
  # L^9 + 2L^8 - L^6 - 2L^5 + L^3 [Example (4) in Section 3.2, loc. cit.]
  @test string(nullcone_motive(jordan_quiver(2), [3])) ==
    "L^9 + 2*L^8 - L^6 - 2*L^5 + L^3"

  # the defect measures the failure of equidimensionality: it vanishes for cofree
  # settings, and for the settings 2 <=> k it decreases to zero as k grows to 4
  @test defect(jordan_quiver(1), [4]) == 0
  @test defect(jordan_quiver(2), [2]) == 0
  @test defect(jordan_quiver(2), [3]) == 1
  @test [defect(Quiver("1--2, 2-1"), [2, k]) for k in 2:4] == [2, 1, 0]
  @test defect(cyclic_quiver(3), [2, 3, 4]) == 0

  # Popov: cofree iff coregular with vanishing defect, as an independent check of
  # the classification-based is_cofree against the independently computed defect
  for a in 0:2, b in 0:2, l1 in 0:1, l2 in 0:1, d1 in 1:2, d2 in 1:2
    a + b + l1 + l2 <= 3 || continue
    Q = Quiver([l1 a; b l2])
    d = [d1, d2]
    @test is_cofree(Q, d) == (is_coregular(Q, d) && defect(Q, d) == 0)
  end
end;
QuiverTools.empty_all_caches!()

@testset "projections to walls             " begin
  # the projection of the 6-subspace quiver moduli with d = (1^5, 2; 3) from the
  # canonical chamber to the wall (1^5, 2; -3), whose target is Gr(2, 4): birational,
  # with fibres P^1 over five surfaces and P^1 x P^1 over ten points, hence semismall
  # but not flat
  S = subspace_quiver(6)
  d = [1, 1, 1, 1, 1, 2, 3]
  theta = [3, 3, 3, 3, 3, 3, -7]
  thetabar = [1, 1, 1, 1, 1, 2, -3]
  @test fibre_dimension(S, d, theta, thetabar, Dict(d => [1])) == 0
  e1, e3 = [1, 0, 0, 0, 0, 1, 1], [0, 1, 1, 1, 1, 1, 2]
  @test fibre_dimension(S, d, theta, thetabar, Dict(e1 => [1], e3 => [1])) == 1
  e2, f2 = [0, 1, 0, 0, 0, 1, 1], [0, 0, 1, 1, 1, 0, 1]
  @test fibre_dimension(S, d, theta, thetabar, Dict(e1 => [1], e2 => [1], f2 => [1])) == 2
  @test !is_flat(S, d, theta, thetabar)
  @test is_semismall(S, d, theta, thetabar)

  # the projection to the other wall is a semismall resolution of the ten isolated
  # singularities, with fibres P^2; the target is singular so flatness would need more
  # than miracle flatness
  @test is_semismall(S, d, theta, [2, 2, 2, 2, 2, 1, -4])
  @test_throws ArgumentError is_flat(S, d, theta, [2, 2, 2, 2, 2, 1, -4])

  # the identity projection is flat
  @test is_flat(kronecker_quiver(3), [2, 3], [3, -2], [3, -2])

  # the source stability parameter must be King-normalized
  @test_throws ArgumentError fibre_dimension(
    S, d, [1, 1, 1, 1, 1, 1, -7], thetabar, Dict(d => [1])
  )

  # The target must also be King-normalized, and the source must be a chamber point.
  @test_throws ArgumentError fibre_dimension(
    S, d, theta, thetabar .+ [1, 0, 0, 0, 0, 0, 0], Dict(d => [1])
  )
  @test_throws ArgumentError fibre_dimension(S, d, thetabar, thetabar, Dict(d => [1]))

  # Parameters in different VGIT chambers do not define a projection to a wall.
  T = three_vertex_quiver(2, 3, 4)
  @test_throws ArgumentError fibre_dimension(
    T, [1, 2, 2], [4, -1, -1], [4, 1, -3], Dict([1, 2, 2] => [1])
  )

  # A custom denominator belongs to the original quiver only. This closure accepts
  # seven-component vectors and would fail if reused on the two-vertex local quiver.
  weights = collect(1:7)
  weighted_denom = e -> sum(weights .* e)
  @test fibre_dimension(
    S, d, theta, thetabar, Dict(e1 => [1], e3 => [1]); denom=weighted_denom
  ) == 1
  @test_throws ArgumentError fibre_dimension(
    S, d, theta, thetabar, Dict(d => [1]); denom=_e -> 0
  )
end;
QuiverTools.empty_all_caches!()

@testset "Domokos reductions               " begin
  # M. Domokos, Quiver moduli spaces of a given dimension, https://doi.org/10.4171/JCA/97
  # τ and σ reductions preserve the moduli space up to isomorphism (Theorem 2.5), which we
  # witness through the Hodge diamond. The stability sign is QuiverTools' (opposite King),
  # so the Section 9 weight θ_paper = (-2, 1, 1) is passed as (2, -1, -1).
  function diag(M)
    h = hodge_diamond(M)
    return [h[i, i] for i in 1:size(h, 1)]
  end

  # τ_u at a large middle vertex: 1 → 2 ⇉⇉⇉ 3 reduces to the 3-Kronecker quiver.
  # With d = [1, 3, 1] the moduli space is ℙ²; with d = [2, 9, 3] it is the 6-fold M(K₃, [2,3]).
  Q = Quiver("1-2,2---3")
  @test is_large(Q, [1, 3, 1], 2)
  @test tau_reduction(Q, [1, 3, 1], 2) == (Quiver([0 3; 0 0]), [1, 1])
  @test tau_reduction(Q, [1, 3, 1], [3, 1, -6], 2) ==
    (Quiver([0 3; 0 0]), [1, 1], [3, -3])
  @test diag(QuiverModuliSpace(Q, [1, 3, 1], [3, 1, -6])) == [1, 1, 1]              # ℙ²
  @test diag(QuiverModuliSpace(Q, [2, 9, 3], [9, 1, -9])) ==
    diag(tau_reduction(QuiverModuliSpace(Q, [2, 9, 3], [9, 1, -9]), 2)) ==
    [1, 1, 3, 3, 3, 1, 1]                                                           # the 6-fold

  # Section 9 (Proposition 9.1): a τσ-minimal pair with no large vertex.
  Q9, d9, th9 = Quiver("1--3,1-2,3-2"), [2, 1, 3], [2, -1, -1]
  @test !any(v -> is_large(Q9, d9, v), 1:3)
  @test is_small_source(Q9, d9, 1) && is_small_sink(Q9, d9, 2)
  @test tau_sigma_reduce(Q9, d9) == (Q9, d9)
  @test tau_sigma_reduce(Q9, d9, th9) == (Q9, d9, th9)          # greedy reduce is a fixpoint
  # σ preserves the moduli space (dimension is a cheap invariant here)
  M9 = QuiverModuliSpace(Q9, d9, th9)
  @test dimension(sigma_reduction(M9, 1)) == dimension(M9) == 4

  # σ_u is an involution
  @test sigma_reduction(sigma_reduction(Q9, d9, 1)..., 1) == (Q9, d9)
  @test sigma_reduction(sigma_reduction(Q9, d9, th9, 1)..., 1) == (Q9, d9, th9)

  # the 2-cycle with d = [1, 1]: not minimal (large vertex), reduces to a single loop
  @test is_taus_minimal(Quiver("1-2,2-1"), [1, 1]) == false
  @test tau_sigma_reduce(Quiver("1-2,2-1"), [1, 1]) == (Quiver([1;;]), [1])
  @test tau_sigma_reduce(Quiver("1-2,2-1"), [1, 1], [1, -1]) == (Quiver([1;;]), [1], [0])
  @test is_taus_minimal(Quiver([1;;]), [1]) == true

  # A denominator closure belongs to its original vertex set. The transformed King
  # weight is normalized, so reduced moduli objects use the standard denominator.
  weighted_denom = e -> [1, 2, 3]' * e
  @test tau_reduction(
    QuiverModuliSpace(Q, [1, 3, 1], [3, 1, -6], "semistable", weighted_denom), 2
  ).denom === sum
  @test sigma_reduction(
    QuiverModuliSpace(Q9, d9, th9, "semistable", weighted_denom), 1
  ).denom === sum

  # error paths: the vertex must be large resp. small, theta must be King-normalized,
  # the source must be semistable, and the dimension vector must be sincere
  @test_throws ArgumentError tau_reduction(Q, [1, 3, 1], [3, 1, -6], 1)
  @test_throws ArgumentError tau_reduction(Q, [1, 3, 1], [3, 1, -5], 2)
  @test_throws ArgumentError tau_reduction(Quiver("1-2,2-1"), [2, 1], [1, -2], 1)
  @test_throws ArgumentError sigma_reduction(Q9, d9, th9, 3)
  @test_throws ArgumentError sigma_reduction(Q9, d9, [2, -1, 0], 1)
  @test_throws ArgumentError tau_sigma_reduce(Quiver([1;;]), [1], [1])
  @test_throws AssertionError is_large(Q, [1, 0, 1], 2)
  @test_throws AssertionError is_taus_minimal(Quiver([1;;]), [0])
  @test_throws ArgumentError is_taus_minimal(Quiver([1;;]), [1]; max_states=0)

  # The sigma-orbit of the Section 9 pair is infinite, so a bounded search cannot prove
  # its minimality. Hitting the cap reports an indeterminate result rather than `true`,
  # and a still larger traversal reports integer overflow rather than mistaking wrapped
  # dimension entries for a smaller pair.
  @test_throws ArgumentError is_taus_minimal(Q9, d9; max_states=2)
  @test_throws OverflowError is_taus_minimal(Q9, d9; max_states=200)

  # With this pair the first reflection does not lower the size, but it exposes a later
  # tau reduction. A cap of one therefore must not produce a false proof of minimality.
  Q4 = Quiver([0 0 1 0; 0 0 2 2; 0 1 0 1; 0 2 1 0])
  d4 = [3, 2, 6, 1]
  @test_throws ArgumentError is_taus_minimal(Q4, d4; max_states=1)
  @test is_taus_minimal(Q4, d4; max_states=100) == false
end;
QuiverTools.empty_all_caches!()

# The number of `jldoctest` blocks in the docstrings of `mod`. `doctest` runs them all
# as a single test, so the number is shown in the name of the testset instead.
function count_doctests(mod::Module)
  return sum(values(Base.Docs.meta(mod)); init=0) do multidoc
    sum(values(multidoc.docs); init=0) do docstr
      count("```jldoctest", join(docstr.text))
    end
  end
end

@testset "$(rpad("Doctests ($(count_doctests(QuiverTools)) blocks)", 33))" begin
  # `import Oscar` (not `using`) loads Oscar so the walls-and-chambers / VGIT extension
  # activates and its doctests resolve `Oscar.*`, without pulling Oscar's exports into
  # scope (which would clash with QuiverTools names such as `index`, `todd_class`, ...).
  DocMeta.setdocmeta!(
    QuiverTools, :DocTestSetup, :(using QuiverTools; import Oscar, GraphViz)
  )
  # show only Documenter's warnings and errors, which include doctest failures
  with_logger(ConsoleLogger(stderr, Logging.Warn)) do
    doctest(QuiverTools; manual=true, testset="Doctests")
  end
end;
