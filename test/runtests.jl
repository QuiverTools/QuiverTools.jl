using Test, QuiverTools, Documenter
using Pkg;
Pkg.activate(@__DIR__)

@info "Almost all the tests are in the documentation."

DocMeta.setdocmeta!(QuiverTools, :DocTestSetup, :(using QuiverTools))
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
  open("3vertexquiver-3-5-7-canonical.txt", "r") do file
    expected = readline(file)
  end

  @test string(all_hn_types(Q, d, theta; ordered=true)) == expected
end;

@testset "covering quiver" begin
  # Argument validation
  @test_throws ArgumentError compatible_dimension_vectors(kronecker_quiver(2), [1, 1, 1])
  @test_throws ArgumentError compatible_dimension_vectors(kronecker_quiver(2), [-1, 1])

  # Zero dimension vector: one (empty) shift class
  @test compatible_dimension_vectors(kronecker_quiver(2), [0, 0]) == [CoveringDimVector()]

  @testset "Kronecker counts" begin
    # m-Kronecker with d = (1, 1) gives P^{m-1} with m torus-fixed points.
    # Cited in Section 3 of arXiv:2002.12049 as the prototypical example.
    for m in 1:5
      Q = kronecker_quiver(m)
      @test length(compatible_dimension_vectors(Q, [1, 1])) == m
    end

    # 3-Kronecker, d = (2, 3): 55 shift classes
    @test length(compatible_dimension_vectors(kronecker_quiver(3), [2, 3])) == 55

    # 3-Kronecker, d = (1, 2): a smaller non-trivial case
    @test length(compatible_dimension_vectors(kronecker_quiver(3), [1, 2])) == 6
  end

  @testset "covering_euler_form = euler_form on induced subquiver" begin
    # For any finitely-supported beta, the covering Euler form coincides with
    # the underlying Euler form on the induced finite subquiver.
    for (Q, d) in [
      (kronecker_quiver(2), [1, 1]),
      (kronecker_quiver(3), [1, 1]),
      (kronecker_quiver(3), [2, 3]),
      (subspace_quiver(3), [1, 1, 1, 2]),
    ]
      for beta in compatible_dimension_vectors(Q, d)
        sub_Q, sub_d, _ = extract_finite_subquiver(Q, beta)
        @test covering_euler_form(Q, beta, beta) ==
          euler_form(sub_Q, sub_d, sub_d)
      end
    end
  end

  @testset "shift invariance" begin
    # The Z^{Q_1}-shift action preserves the covering Euler form.
    Q = kronecker_quiver(3)
    beta = CoveringDimVector((1, [0, 0, 0]) => 1, (2, [1, 0, 0]) => 1)
    gamma = CoveringDimVector((1, [0, 0, 0]) => 1, (2, [0, 1, 0]) => 1)
    for chi in ([0, 0, 0], [1, 0, 0], [2, -1, 3], [-5, 4, 7])
      @test covering_euler_form(Q, beta, gamma) ==
        covering_euler_form(Q, shift_beta(beta, chi), shift_beta(gamma, chi))
    end

    # Shift is an involutive group action: s_chi composed with s_{-chi} is identity.
    for chi in ([0, 0, 0], [3, -2, 1])
      @test shift_beta(shift_beta(beta, chi), .-chi) == beta
    end
  end

  @testset "Boos--Franzen Thm 6.1 sum identity (algebraic)" begin
    # Algebraic sum of weight_space_dimension over all chi recovers
    # 1 - <d,d>_Q for every compatible beta:
    # [[Theorem 6.1, arXiv:2002.12049]] gives
    # dim (T_[M])_chi = delta(chi,0) - <beta, s_{-chi} beta>_{Q(w)}, and
    # sum over chi of <beta, s_{-chi} beta>_{Q(w)} collapses to
    # sum_i d_i^2 - sum_a d_{s(a)} d_{t(a)} = <d,d>_Q.
    for (Q, d) in [
      (kronecker_quiver(2), [1, 1]),
      (kronecker_quiver(3), [1, 1]),
      (kronecker_quiver(3), [1, 2]),
      (kronecker_quiver(3), [2, 3]),
      (subspace_quiver(3), [1, 1, 1, 2]),
    ]
      expected = 1 - euler_form(Q, d, d)
      for beta in compatible_dimension_vectors(Q, d)
        candidates = QuiverTools._weight_candidates(Q, beta)
        total = sum(weight_space_dimension(Q, beta, chi) for chi in candidates)
        @test total == expected
      end
    end
  end

  @testset "real-root fixed points are smooth" begin
    # At a real-root beta (covering Euler form 1), the fixed-point component
    # is isolated and all tangent weight space dimensions are nonnegative.
    # In particular, the sum of nonzero weight dimensions agrees with
    # 1 - <d,d>_Q.
    for (Q, d) in [
      (kronecker_quiver(2), [1, 1]),
      (kronecker_quiver(3), [1, 1]),
      (kronecker_quiver(4), [1, 1]),
    ]
      expected = 1 - euler_form(Q, d, d)
      for beta in compatible_dimension_vectors(Q, d)
        covering_euler_form(Q, beta, beta) == 1 || continue
        weights = nonzero_weights(Q, beta)
        @test all(dim > 0 for (_, dim) in weights)
        @test sum(dim for (_, dim) in weights; init=0) == expected
      end
    end
  end

  @testset "Kronecker P^1 fixed point partners" begin
    # The two torus-fixed points of P^1 = M^theta(K_2, (1,1)) have nonzero
    # tangent weights that are negatives of each other.
    Q = kronecker_quiver(2)
    betas = compatible_dimension_vectors(Q, [1, 1])
    weights = [nonzero_weights(Q, beta) for beta in betas]
    @test length(weights) == 2
    @test all(length(w) == 1 for w in weights)
    chi1 = first(weights[1])[1]
    chi2 = first(weights[2])[1]
    @test chi1 == .-chi2
  end

  @testset "3-Kronecker P^2 valency" begin
    # P^2 = M^theta(K_3, (1,1)) has 3 torus-fixed points; at each, the tangent
    # space is 2-dimensional, decomposing into 2 weight spaces of dimension 1
    # whose characters are pairwise Q-linearly independent.
    Q = kronecker_quiver(3);
    d = [1, 1]
    betas = compatible_dimension_vectors(Q, d)
    @test length(betas) == 3
    for beta in betas
      weights = nonzero_weights(Q, beta)
      @test length(weights) == 2
      @test all(dim == 1 for (_, dim) in weights)
      chi1, chi2 = weights[1][1], weights[2][1]
      # Q-linear independence: not a rational multiple.
      @test any(
        chi1[i] * chi2[j] != chi1[j] * chi2[i]
        for i in eachindex(chi1), j in eachindex(chi1)
      )
    end
  end

  @testset "canonical form is idempotent" begin
    # compatible_dimension_vectors returns canonical representatives; if we
    # extract one and renormalise it (via a no-op shift), we get the same thing.
    Q = kronecker_quiver(3)
    for beta in compatible_dimension_vectors(Q, [2, 3])
      @test shift_beta(beta, zeros(Int, n_arrows(Q))) == beta
    end
  end
end;

# @testset "Testing weight handling" begin

#     U = Bundle([1,2],2)
#     V = Bundle([3,3,3],3)

#     @test U ⊕ V == Bundle([1, 2, 3, 3, 3], 5)
#     @test U ⊗ V == Bundle([4, 4, 4, 5, 5, 5], 6)
#     @test U ⊠ V == Bundle([ [1, 3], [1, 3], [1, 3],
#                             [2, 3], [2, 3], [2, 3]], 6)

#     @test U - V == Bundle([1, 2, -3, -3, -3], 5)

#     @test wedge(U,0) == Bundle([0], 1)
#     @test wedge(U,1) == Bundle([1, 2], 2)
#     @test wedge(U,2) == Bundle([3], 1)
#     @test wedge(U,3) == Bundle(Int64[], 0)
#     @test wedge(U ⊠ V,2) == Bundle([[2, 6], [2, 6], [3, 6],
#                                     [3, 6], [3, 6], [2, 6],
#                                     [3, 6], [3, 6], [3, 6],
#                                     [3, 6], [3, 6], [3, 6],
#                                     [4, 6], [4, 6], [4, 6]], 15)

#     @test wedge(U ⊠ V,3) == Bundle([[3, 9], [4, 9], [4, 9],
#                                     [4, 9], [4, 9], [4, 9],
#                                     [4, 9], [5, 9], [5, 9],
#                                     [5, 9], [4, 9], [4, 9],
#                                     [4, 9], [5, 9], [5, 9],
#                                     [5, 9], [5, 9], [5, 9],
#                                     [5, 9], [6, 9]], 20)

#     @test wedge(U ⊠ V,4) == Bundle([[5, 12], [5, 12], [5, 12],
#                                     [6, 12], [6, 12], [6, 12],
#                                     [6, 12], [6, 12], [6, 12],
#                                     [7, 12], [6, 12], [6, 12],
#                                     [6, 12], [7, 12], [7, 12]], 15)

#     @test wedge(U ⊠ V,5) == Bundle([[7, 15], [7, 15], [7, 15],
#                                     [8, 15], [8, 15], [8, 15]], 6)

#     @test wedge(U ⊠ V,6) == Bundle([[9, 18]], 1)

#     W = Bundle([1, 2, 3], 3);

#     @test symm(W,0) == Bundle([0], 1)
#     @test symm(W,1) == Bundle([1, 2, 3], 3)
#     @test symm(W,2) == Bundle([2, 3, 4, 4, 5, 6], 6)
#     @test symm(W,3) == Bundle([3, 4, 5, 5, 6, 7, 6, 7, 8, 9], 10)
#     @test symm(W,4) == Bundle([4, 5, 6, 6, 7, 8, 7, 8, 9, 10, 8, 9, 10, 11, 12], 15)
# end;
