@testset "covering quiver" begin
  @testset "input validation" begin
    Q = kronecker_quiver(2)
    beta = CoveringDimVector((1, [0, 0]) => 1)

    @test_throws ArgumentError compatible_dimension_vectors(Q, [1, 1, 1])
    @test_throws ArgumentError compatible_dimension_vectors(Q, [-1, 1])
    @test_throws ArgumentError extract_finite_subquiver(Q, beta, [1])
    @test_throws ArgumentError weight_space_dimension(Q, beta, [0])
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
    @test issorted(first.(nonzero_weights(kronecker_quiver(3), beta)))
  end

  @testset "quivers without arrows" begin
    Q = Quiver(zeros(Int, 1, 1))
    beta = only(compatible_dimension_vectors(Q, [2]))

    @test beta == CoveringDimVector((1, Int[]) => 2)
    @test covering_euler_form(Q, beta, beta) == 4
    @test isempty(nonzero_weights(Q, beta))
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
  end

  @testset "tangent weights" begin
    examples = [
      (kronecker_quiver(2), [1, 1]),
      (kronecker_quiver(3), [1, 1]),
      (kronecker_quiver(3), [1, 2]),
      (kronecker_quiver(3), [2, 3]),
      (subspace_quiver(3), [1, 1, 1, 2]),
    ]

    for (Q, d) in examples, beta in compatible_dimension_vectors(Q, d)
      total = sum(
        weight_space_dimension(Q, beta, chi)
        for chi in QuiverTools._weight_candidates(Q, beta)
      )
      @test total == 1 - euler_form(Q, d, d)
    end

    for m in 2:4
      Q = kronecker_quiver(m)
      expected = 1 - euler_form(Q, [1, 1], [1, 1])
      for beta in compatible_dimension_vectors(Q, [1, 1])
        weights = nonzero_weights(Q, beta)
        @test all(last(weight) > 0 for weight in weights)
        @test sum(last, weights; init=0) == expected
      end
    end

    Q = kronecker_quiver(2)
    weights = nonzero_weights.(Ref(Q), compatible_dimension_vectors(Q, [1, 1]))
    @test length(weights) == 2
    @test all(length(weight) == 1 for weight in weights)
    @test only(first.(weights[1])) == .-only(first.(weights[2]))
  end

  @testset "fixed components" begin
    M = QuiverModuliSpace(kronecker_quiver(2), [1, 1], [1, -1])
    components = torus_fixed_components(M)
    @test length(components) == 2
    @test all(component.moduli.condition == "stable" for component in components)
    @test all(dimension(component.moduli) == 0 for component in components)
    @test all(
      sum(values(component.beta); init=0) == sum(M.d; init=0)
      for component in components
    )

    Q = Quiver([0 1; 0 0])
    M = QuiverModuliSpace(Q, [1, 2], [2, -1])
    @test isempty(torus_fixed_components(M))

    M = QuiverModuliSpace(kronecker_quiver(3), [2, 2], [1, -1])
    @test_throws ArgumentError torus_fixed_components(M)

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
  end
end
