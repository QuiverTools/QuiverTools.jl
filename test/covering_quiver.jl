@testset "covering quiver" begin
  @testset "input validation" begin
    Q = kronecker_quiver(2)
    beta = CoveringDimVector((1, [0, 0]) => 1)
    M = QuiverModuliSpace(Q, [1, 1], [1, -1])
    stable_beta = CoveringDimVector(
      (1, [0, 0]) => 1,
      (2, [1, 0]) => 1,
    )

    @test_throws ArgumentError compatible_dimension_vectors(Q, [1, 1, 1])
    @test_throws ArgumentError compatible_dimension_vectors(Q, [-1, 1])
    @test_throws ArgumentError extract_finite_subquiver(Q, beta, [1])
    @test_throws ArgumentError weight_space_dimension(M, stable_beta, [0])
    @test_throws ArgumentError weight_space_dimension(M, beta, [0, 0])

    invalid_vertex = CoveringDimVector((3, [0, 0]) => 1)
    @test_throws ArgumentError extract_finite_subquiver(Q, invalid_vertex)
    @test_throws ArgumentError covering_euler_form(Q, invalid_vertex, beta)

    nonpositive_vertex = CoveringDimVector((0, [0, 0]) => 1)
    @test_throws ArgumentError extract_finite_subquiver(Q, nonpositive_vertex)

    wrong_coordinate = CoveringDimVector((1, [0]) => 1)
    @test_throws DimensionMismatch extract_finite_subquiver(Q, wrong_coordinate)
    @test_throws DimensionMismatch shift_beta(beta, [0])

    zero_multiplicity = CoveringDimVector((1, [0, 0]) => 0)
    negative_multiplicity = CoveringDimVector((1, [0, 0]) => -1)
    @test_throws ArgumentError extract_finite_subquiver(Q, zero_multiplicity)
    @test_throws ArgumentError shift_beta(negative_multiplicity, [0, 0])

    gamma = CoveringDimVector((2, [0]) => 1)
    @test_throws DimensionMismatch covering_euler_form(Q, beta, gamma)

    empty_lift = QuiverModuliSpace(Quiver([0 1; 0 0]), [1, 2], [2, -1], "stable")
    empty_beta = CoveringDimVector((1, [0]) => 1, (2, [1]) => 2)
    @test_throws ArgumentError weight_space_dimension(empty_lift, empty_beta, [0])
  end

  @test compatible_dimension_vectors(kronecker_quiver(2), [0, 0]) ==
    [CoveringDimVector()]

  @testset "Kronecker counts" begin
    for m in 1:5
      @test length(compatible_dimension_vectors(kronecker_quiver(m), [1, 1])) == m
    end

    @test length(compatible_dimension_vectors(kronecker_quiver(3), [1, 2])) == 6
    @test length(compatible_dimension_vectors(kronecker_quiver(3), [2, 3])) == 55
  end

  @testset "deterministic representatives" begin
    betas = compatible_dimension_vectors(kronecker_quiver(3), [2, 3])
    canonical_betas = QuiverTools._canonicalize.(betas)
    @test issorted(canonical_betas)

    beta = CoveringDimVector(
      (1, [0, 0, 0]) => 1,
      (2, [1, 0, 0]) => 1,
    )
    M = QuiverModuliSpace(kronecker_quiver(3), [1, 1], [1, -1])
    @test issorted(first.(tangent_weight_multiplicities(M, beta)))
  end

  @testset "quivers without arrows" begin
    Q = Quiver(zeros(Int, 1, 1))
    beta = only(compatible_dimension_vectors(Q, [2]))

    @test beta == CoveringDimVector((1, Int[]) => 2)
    @test covering_euler_form(Q, beta, beta) == 4
    @test shift_beta(beta, Int[]) == beta
  end

  @testset "finite support calculations" begin
    examples = [
      (kronecker_quiver(2), [1, 1]),
      (kronecker_quiver(3), [1, 1]),
      (kronecker_quiver(3), [2, 3]),
      (subspace_quiver(3), [1, 1, 1, 2]),
    ]

    for (Q, d) in examples, beta in compatible_dimension_vectors(Q, d)
      sub_Q, sub_d, _ = extract_finite_subquiver(Q, beta)
      @test covering_euler_form(Q, beta, beta) == euler_form(sub_Q, sub_d, sub_d)
    end

    Q = kronecker_quiver(3)
    beta = CoveringDimVector(
      (1, [0, 0, 0]) => 1,
      (2, [1, 0, 0]) => 1,
    )
    gamma = CoveringDimVector(
      (1, [0, 0, 0]) => 1,
      (2, [0, 1, 0]) => 1,
    )
    for chi in ([0, 0, 0], [1, 0, 0], [2, -1, 3], [-5, 4, 7])
      @test covering_euler_form(Q, beta, gamma) == covering_euler_form(
        Q,
        shift_beta(beta, chi),
        shift_beta(gamma, chi),
      )
    end
    @test shift_beta(shift_beta(beta, [3, -2, 1]), [-3, 2, -1]) == beta

    _, sub_d, sub_theta, vertex_map = extract_finite_subquiver(Q, beta, [1, -1])
    vertices = sort(collect(keys(vertex_map)))
    @test sub_d == [beta[vertex] for vertex in vertices]
    @test sub_theta == [vertex[1] == 1 ? 1 : -1 for vertex in vertices]

    pushed_beta = [1, 1]
    pushed_gamma = [1, 1]
    shifted_sum = sum(
      covering_euler_form(Q, beta, shift_beta(gamma, [-i, -j, -k])) for
      i in -2:2 for j in -2:2 for k in -2:2
    )
    @test shifted_sum == euler_form(Q, pushed_beta, pushed_gamma)
  end

  @testset "tangent weights" begin
    examples = QuiverModuliSpace[
      QuiverModuliSpace(kronecker_quiver(2), [1, 1], [1, -1]),
      QuiverModuliSpace(kronecker_quiver(3), [1, 1], [1, -1]),
      QuiverModuliSpace(kronecker_quiver(3), [2, 3]),
      QuiverModuliSpace(subspace_quiver(3), [1, 1, 1, 2]),
    ]

    for M in examples
      components = torus_fixed_components(M)
      for component in components
        weights = tangent_weight_multiplicities(M, component.beta)
        @test all(last(weight) > 0 for weight in weights)
        @test sum(last, weights; init=0) == 1 - euler_form(M.Q, M.d, M.d)
      end
    end

    for m in 2:4
      Q = kronecker_quiver(m)
      M = QuiverModuliSpace(Q, [1, 1], [1, -1])
      expected = 1 - euler_form(Q, [1, 1], [1, 1])
      for component in torus_fixed_components(M)
        weights = tangent_weight_multiplicities(M, component.beta)
        @test all(last(weight) > 0 for weight in weights)
        @test sum(last, weights; init=0) == expected
      end
    end

    Q = kronecker_quiver(2)
    M = QuiverModuliSpace(Q, [1, 1], [1, -1])
    weights = [
      tangent_weight_multiplicities(M, component.beta) for
      component in torus_fixed_components(M)
    ]
    @test length(weights) == 2
    @test all(length(weight) == 1 for weight in weights)
    @test only(first.(weights[1])) == .-only(first.(weights[2]))

    Q = kronecker_quiver(3)
    M = QuiverModuliSpace(Q, [3, 4], [4, -3], "stable")
    positive_dimensional_beta = CoveringDimVector(
      (1, [0, 0, 0]) => 1,
      (1, [0, 1, -1]) => 1,
      (1, [1, 0, -1]) => 1,
      (2, [0, 0, 1]) => 1,
      (2, [0, 1, 0]) => 1,
      (2, [1, 0, 0]) => 1,
      (2, [1, 1, -1]) => 1,
    )
    @test weight_space_dimension(M, positive_dimensional_beta, zeros(Int, 3)) == 1
    weights = tangent_weight_multiplicities(M, positive_dimensional_beta)
    @test ([0, 0, 0], 1) in weights
    @test sum(last, weights; init=0) == 1 - euler_form(Q, M.d, M.d)

    Q = jordan_quiver()
    M = QuiverModuliSpace(Q, [1], [0])
    component = only(torus_fixed_components(M))
    @test tangent_weight_multiplicities(M, component.beta) == [([1], 1)]
  end

  @testset "fixed components" begin
    M = QuiverModuliSpace(kronecker_quiver(2), [1, 1], [1, -1])
    components = torus_fixed_components(M)
    @test length(components) == 2
    @test all(component.moduli.condition == "stable" for component in components)
    @test all(dimension(component.moduli) == 0 for component in components)
    @test all(
      [
        sum(
          multiplicity for ((v, _), multiplicity) in component.beta if v == i;
          init=0,
        ) for i in eachindex(M.d)
      ] == M.d
      for component in components
    )

    Q = Quiver([0 1; 0 0])
    M = QuiverModuliSpace(Q, [1, 2], [2, -1])
    @test isempty(torus_fixed_components(M))

    M = QuiverModuliSpace(kronecker_quiver(3), [2, 3])
    @test length(torus_fixed_components(M)) == 13

    M = QuiverModuliSpace(kronecker_quiver(3), [2, 2], [1, -1])
    @test_throws ArgumentError torus_fixed_components(M)
    @test_throws ArgumentError tangent_weight_multiplicities(
      M,
      first(compatible_dimension_vectors(M.Q, M.d)),
    )

    nonnormalized = QuiverModuliSpace(
      Quiver(zeros(Int, 1, 1)),
      [2],
      [1],
      "semistable",
    )
    @test_throws ArgumentError torus_fixed_components(nonnormalized)

    stable_locus = QuiverModuliSpace(
      kronecker_quiver(3),
      [2, 2],
      [1, -1],
      "stable",
    )
    @test torus_fixed_components(stable_locus) isa Vector

    denom = d -> 2 * d[1] + d[2]
    custom_denom_locus = QuiverModuliSpace(
      kronecker_quiver(3),
      [1, 2],
      [2, -1],
      "stable",
      denom,
    )
    components = torus_fixed_components(custom_denom_locus)
    @test !isempty(components)
    @test all(
      component.moduli.denom(component.moduli.d) == denom(custom_denom_locus.d)
      for component in components
    )
    for component in components
      _, _, vertex_map = extract_finite_subquiver(custom_denom_locus.Q, component.beta)
      for ((v, _), sub_v) in vertex_map
        sub_e = zeros(Int, length(vertex_map))
        sub_e[sub_v] = 1
        projected_e = zeros(Int, n_vertices(custom_denom_locus.Q))
        projected_e[v] = 1
        @test component.moduli.denom(sub_e) == denom(projected_e)
      end
    end
  end
end
