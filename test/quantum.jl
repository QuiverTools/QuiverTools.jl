@testset "Quasimap quantum cohomology" begin
  S = QuiverTools.Singular
  N = S.Nemo
  normal(f) = QuiverTools.__simplify(f)
  finite_rank(R) = length(QuiverTools.gens(S.kbase(QuiverTools.quotient_ideal(R))))
  classical(A) = specialize_quantum_parameters(A, zeros(Int, length(novikov_variables(A))))

  # Independent reference for the reduced-word implementation: explicitly sum
  # over the Weyl group and divide by its discriminant.
  function explicit_alternant(f, d, roots)
    R = parent(f)
    offsets = cumsum([0; d])
    blocks = [collect((offsets[i] + 1):offsets[i + 1]) for i in eachindex(d)]
    delta = prod(
      roots[j] - roots[i]
      for block in blocks for i in block for j in block if i < j; init=R(1))
    alternating = R(0)
    for w in Iterators.product([QuiverTools.permutations(block) for block in blocks]...)
      permutation = [vcat(w...); collect((length(roots) + 1):length(QuiverTools.gens(R)))]
      weyl_sign = prod(sign(S.perm(sortperm(block))) for block in w; init=1)
      context = S.MPolyBuildCtx(R)
      for (c, e) in zip(S.coefficients(f), S.exponent_vectors(f))
        S.push_term!(context, weyl_sign * c, e[permutation])
      end
      alternating += S.finish(context)
    end
    return S.divexact(alternating, delta)
  end

  @testset "Reduced-word divided alternants" begin
    for d in ([1], [0, 2, 1], [3], [2, 3])
      r = sum(d)
      R, vars = S.polynomial_ring(S.QQ, ["z$i" for i in 1:r] ∪ ["q1", "q2"])
      roots, q = vars[1:r], vars[(r + 1):end]
      differences, basis = QuiverTools.__qm_weyl_data(d, R, roots)
      @test length(differences) == sum(n * (n - 1) ÷ 2 for n in d)
      f =
        sum(i * roots[i]^(i + 1) for i in 1:r; init=R(0)) +
        q[1] * prod(roots[i]^(i - 1) for i in 1:r; init=R(1)) + q[2]^2
      for b in basis
        @test QuiverTools.__qm_divided_alternant(b * f, differences) ==
          explicit_alternant(b * f, d, roots)
      end
      @test iszero(QuiverTools.__qm_divided_alternant(R(0), differences))
    end
  end

  @testset "Projective spaces and nilpotents" begin
    for n in 2:5
      M = QuiverModuliSpace(kronecker_quiver(n), [1, 1])
      A = quasimap_quantum_cohomology(M)
      x, q = quantum_tautological_classes(A), only(novikov_variables(A))
      h = x[2][1] - x[1][1]
      @test iszero(normal(h^n - q))
      @test iszero(quantum_normal_form(A, h^n - q))
      @test gw_comparison(A) == :equal
      @test finite_rank(classical(A)) == n
      matrix, basis, R = quantum_multiplication_matrix(A, h, [1])
      @test matrix^n == one(matrix)
      @test length(basis) == n
      nilpotent, _, _ = quantum_multiplication_matrix(A, h, [0])
      @test iszero(nilpotent^n)
      @test !iszero(nilpotent^(n - 1))
      @test quasimap_presentation(A) === A.ring
      @test !isdefined(M.chow, :ring)
      degrees = novikov_degrees(A)
      degrees[1][1] = 100
      @test novikov_degrees(A)[1][1] != 100
    end
  end

  @testset "Grassmannian quantum Schur relations" begin
    # Gr(2,4): h_3(c1,c2)=0, h_4(c1,c2)=-q. This checks the Weyl sign.
    A = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(4), [1, 2]);
      chi=[1, 0])
    c1, c2 = quantum_tautological_classes(A)[2]
    q = only(novikov_variables(A))
    @test iszero(normal(c1^3 - 2c1*c2))
    @test iszero(normal(c1^4 - 3c1^2*c2 + c2^2 + q))
    @test finite_rank(classical(A)) == 6
    matrix, _, _ = quantum_multiplication_matrix(A, c1, [1])
    @test matrix^6 == 4matrix^2
    T, t = N.polynomial_ring(N.QQ, "t")
    @test N.charpoly(T, matrix) == t^2 * (t^4 - 4)
    # The zero eigenspace has dimension two; a regular-root solver must not discard it.
    @test N.rank(matrix) == 4
  end

  @testset "Several Novikov parameters" begin
    M = QuiverModuliSpace(three_vertex_quiver(2, 2, 0), [1, 1, 1])
    A = quasimap_quantum_cohomology(M; chi=[1, 0, 0])
    x = quantum_tautological_classes(A)
    q, degrees = novikov_variables(A), novikov_degrees(A)
    @test length(q) == 2
    q_for_vertex(i) = only(q[j] for j in eachindex(q)
                                    if degrees[j][i] - degrees[j][1] == 1)
    h1, h2 = x[2][1], x[3][1]
    @test iszero(normal(h1^2 - q_for_vertex(2)))
    @test iszero(normal(h2^2 - q_for_vertex(3)))
    @test finite_rank(classical(A)) == 4
    @test gw_comparison(A) == :equal
    a, _, _ = quantum_multiplication_matrix(A, h1, [1, 1])
    b, _, _ = quantum_multiplication_matrix(A, h2, [1, 1])
    @test a^2 == one(a)
    @test b^2 == one(b)
    @test a*b == b*a
    T, t = N.polynomial_ring(N.QQ, "t")
    @test N.charpoly(T, a + b) == t^2 * (t^2 - 4)
  end

  @testset "Conservative comparison and presentation dependence" begin
    # F_1 has effective anticanonical degree one: the sufficient test is inconclusive.
    M = QuiverModuliSpace(three_vertex_quiver(2, 1, 1), [1, 1, 1])
    A = quasimap_quantum_cohomology(M; chi=[1, 0, 0])
    x = quantum_tautological_classes(A)
    q, degrees = novikov_variables(A), novikov_degrees(A)
    base = only(q[j] for j in eachindex(q) if degrees[j][2] - degrees[j][1] == 1)
    fibre = only(q[j] for j in eachindex(q) if degrees[j][2] - degrees[j][1] == 0)
    a, b = x[2][1], x[3][1]
    @test iszero(normal(a^2 - base*(b-a)))
    @test iszero(normal(b*(b-a) - fibre))
    @test finite_rank(classical(A)) == 4
    @test gw_comparison(A) == :unknown
    # A forced-nonzero arrow gives a point with a nontrivial quasimap degree monoid.
    B = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(1), [1, 1]))
    @test finite_rank(classical(B)) == 1
    @test length(novikov_variables(B)) == 1
    @test gw_comparison(B) == :unknown
  end

  @testset "Effective degree lattice and monoid kernel" begin
    for d in ([2, 3], [0, 2, 3], [3, 4, 5], [1])
      U = QuiverTools.__qm_degree_lattice(d)
      @test U * BigInt.(d) == [1; zeros(BigInt, length(d)-1)]
      @test abs(N.det(N.matrix(N.ZZ, U))) == 1
    end
    R, q = S.polynomial_ring(S.QQ, ["q1", "q2", "q3"])
    relations = QuiverTools.__qm_monoid_relations(R, q,
      [BigInt[1, 0], BigInt[1, 1], BigInt[1, 2]])
    I = S.std(S.Ideal(R, relations))
    @test iszero(S.reduce(q[1]*q[3] - q[2]^2, I))
    @test !iszero(S.reduce(q[1] - q[2], I))
    @test QuiverTools.__qm_factor(BigInt[3, 2], BigInt(3),
      [BigInt[1, 0], BigInt[1, 1]], BigInt[1, 1]) == [1, 2]
    @test isnothing(
      QuiverTools.__qm_factor(BigInt[1, 2], BigInt(1),
        [BigInt[1, 0], BigInt[1, 1]], BigInt[1, 1]),
    )
  end

  @testset "A nonfree Novikov monoid" begin
    Q = Quiver([0 1 1 0; 0 0 1 1; 0 0 0 1; 0 0 0 0])
    M = QuiverModuliSpace(Q, [1, 1, 1, 2], [6, 3, 1, -5])
    A = quasimap_quantum_cohomology(M)
    q = novikov_variables(A)
    @test length(q) == 4
    @test iszero(normal(q[1]*q[4] - q[2]*q[3]))
    @test finite_rank(classical(A)) == sum(betti_numbers(M)) == 2
    @test_throws ArgumentError specialize_quantum_parameters(A, [1, 1, 1, 2])
    @test finite_rank(specialize_quantum_parameters(A, [1, 1, 1, 1])) == 2
    @test gw_comparison(A) == :unknown
  end

  @testset "Cut effectivity agrees with native stability" begin
    for M in (QuiverModuliSpace(three_vertex_quiver(2, 1, 1), [1, 1, 1]),
      QuiverModuliSpace(kronecker_quiver(3), [1, 2]))
      vertices, edges, cuts = QuiverTools.__qm_cut_data(M)
      r = length(vertices)
      theta = collect(M.theta[vertices])
      for S in QuiverTools.subsets(1:r)
        isempty(S) || length(S) == r ||
          begin
            k = [-Int(i in S) for i in 1:r]
            adjacency = zeros(Int, r, r)
            for (i, j, multiplicity) in edges
              k[j] >= k[i] && (adjacency[i, j] = multiplicity)
            end
            @test (k in cuts) ==
              has_stable_representation(Quiver(adjacency), ones(Int, r), theta)
          end
      end
      # Every effective bounded integral potential splits into effective level cuts.
      for tail in Iterators.product(fill(-2:2, r-1)...)
        k = [0; collect(tail)]
        adjacency = zeros(Int, r, r)
        for (i, j, multiplicity) in edges
          k[j] >= k[i] && (adjacency[i, j] = multiplicity)
        end
        if has_stable_representation(Quiver(adjacency), ones(Int, r), theta)
          levels = [[-Int(a <= level) for a in k] for level in minimum(k):(maximum(k) - 1)]
          @test all(c -> c in cuts, levels)
          @test [sum(c[i] for c in levels; init=0) for i in 1:r] == k .- maximum(k)
        end
      end
    end
  end

  @testset "Nonabelian Kronecker moduli and normalization" begin
    M = QuiverModuliSpace(kronecker_quiver(3), [2, 3])
    A = quasimap_quantum_cohomology(M; chi=[-1, 1])
    @test length(novikov_variables(A)) == 1
    @test finite_rank(classical(A)) == sum(betti_numbers(M)) == 13
    @test finite_rank(specialize_quantum_parameters(A, [1])) == 13
    @test gw_comparison(A) == :equal
    @test !isdefined(M.chow, :ring)
    # A different Bezout normalization changes universal bundles, not the algebra rank.
    B = quasimap_quantum_cohomology(M; chi=[2, -1])
    @test finite_rank(classical(B)) == 13
    x = quantum_tautological_classes(A)
    y = quantum_tautological_classes(B)
    hA, hB = -3x[1][1] + 2x[2][1], -3y[1][1] + 2y[2][1]
    # Meng, arXiv:2412.15987, Theorem 5.23. These degree-two generators
    # are classical in index three, so no change of quantum coordinates is needed.
    c, d, q = x[2][2], x[1][2], only(novikov_variables(A))
    published_relations = [
      hA^4 - 11q*hA + 3c^2 - 9c*d - 3d^2,
      hA^2*c - 3q*hA - c*d - 3d^2,
      hA^2*d - 2q*hA - 3d^2,
      9*(hA*c^2 - q*(4hA^2 - 3c - d)) - 14*(hA*c*d - q*(3hA^2 - c - 3d)),
      2*(hA*c*d - q*(3hA^2 - c - 3d)) - 3*(hA*d^2 - q*(2hA^2 - 3d)),
    ]
    for relation in published_relations
      @test iszero(quantum_normal_form(A, relation))
    end
    a, _, _ = quantum_multiplication_matrix(A, hA, [1])
    b, _, _ = quantum_multiplication_matrix(B, hB, [1])
    T, _ = N.polynomial_ring(N.QQ, "t")
    @test N.charpoly(T, a) == N.charpoly(T, b)
    C = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(4), [2, 3]))
    @test finite_rank(classical(C)) == 58
    @test finite_rank(specialize_quantum_parameters(C, [1])) == 58
    @test gw_comparison(C) == :equal
  end

  @testset "Points, zero support, and validation" begin
    P = quasimap_quantum_cohomology(QuiverModuliSpace(Quiver(zeros(Int, 1, 1)), [1]))
    @test isempty(novikov_variables(P))
    @test finite_rank(classical(P)) == 1
    @test gw_comparison(P) == :equal
    Z = quasimap_quantum_cohomology(
      QuiverModuliSpace(three_vertex_quiver(3, 0, 0),
        [1, 1, 0], [1, -1, 0]),
    )
    @test isempty(quantum_tautological_classes(Z)[3])
    @test finite_rank(classical(Z)) == 3
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(kronecker_quiver(3), [2, 2]))
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(kronecker_quiver(3), [1, 1], [0, 0]))
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(kronecker_quiver(3), [1, 1], [1, 0]))
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(kronecker_quiver(3), [1, 1], [-1, 1]))
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(jordan_quiver(), [1]))
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(Quiver([0 -1; 0 0]), [1, 1], [1, -1]))
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(kronecker_quiver(3), [0, 0]))
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(kronecker_quiver(3), [-1, 1]))
    M = QuiverModuliSpace(kronecker_quiver(3), [1, 1])
    @test_throws ArgumentError quasimap_quantum_cohomology(M; chi=[1, 1])
    @test_throws ArgumentError quasimap_quantum_cohomology(M; chi=[1])
    @test_throws ArgumentError quasimap_quantum_cohomology(M; chi=1)
    @test_throws ArgumentError quasimap_quantum_cohomology(M; chi=[1.0, 0.0])
    @test_throws ArgumentError quasimap_quantum_cohomology(
      QuiverModuliSpace(M.Q, [1, 1], [1, -1], "semistable", x -> sum(x)))
    A = quasimap_quantum_cohomology(M)
    @test_throws ArgumentError specialize_quantum_parameters(A, [])
    @test_throws ArgumentError specialize_quantum_parameters(A, [1.0])
    foreign_ring, _ = S.polynomial_ring(S.QQ, ["z"])
    @test_throws ArgumentError quantum_multiplication_matrix(A, foreign_ring(1), [1])
    @test_throws ArgumentError quantum_normal_form(A, foreign_ring(1))
    @test_throws ArgumentError quantum_normal_form(A, 1)
    @test occursin("Novikov completion", sprint(show, A))
  end
end
