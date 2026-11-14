module ExperimentalChowTests
using Test, QuiverTools
const QT = QuiverTools
const S = QT.Singular

# Independent reference: explicit permutation sum, without any divided differences.
function explicit_alternant(f, d)
  R = parent(f)
  roots = QT.gens(R)
  offsets = cumsum([0; d])
  blocks = [collect((offsets[i] + 1):offsets[i + 1]) for i in eachindex(d)]
  delta = prod(
    roots[j] - roots[i] for block in blocks for i in block for j in block if i < j;
    init=R(1),
  )
  result = R(0)
  for permutation in Iterators.product([QT.permutations(1:r) for r in d]...)
    indices = reduce(vcat, ([offsets[i] + j for j in permutation[i]] for i in eachindex(d)))
    append!(indices, (sum(d) + 1):length(roots))
    sign = prod(S.sign(S.perm(p)) for p in permutation; init=1)
    result += sign * S.AlgebraHomomorphism(R, R, roots[indices])(f)
  end
  S.divexact(result, delta)
end

function transfer(f, R)
  variables = QT.gens(R)
  sum(
    (
      R(c) * prod(variables[i]^e[i] for i in eachindex(variables); init=R(1))
      for (c, e) in zip(S.coefficients(f), S.exponent_vectors(f))
    ); init=R(0))
end

function same_ideal(first_ring, second_ring)
  first_ideal, second_ideal = QT.quotient_ideal(first_ring), QT.quotient_ideal(second_ring)
  first_base, second_base = QT.base_ring(first_ideal), QT.base_ring(second_ideal)
  all(
    iszero(S.reduce(transfer(f, second_base), second_ideal)) for f in QT.gens(first_ideal)
  ) &&
    all(
      iszero(S.reduce(transfer(f, first_base), first_ideal)) for f in QT.gens(second_ideal)
    )
end

@testset "Experimental parabolic Chow backend" begin
  @testset "Operator identities and parabolic module" begin
    for d in ([3], [2, 3], [0, 2, 1])
      r = sum(d)
      R, x = S.polynomial_ring(S.QQ, ["z$i" for i in 1:r] ∪ ["q"])
      f = sum(i * x[i]^(i + 1) for i in 1:r; init=R(0)) + x[end] * x[1]^2
      @test QT.__chow_divided_alternant(f, d) == explicit_alternant(f, d)
    end
    R, x = S.polynomial_ring(S.QQ, ["z$i" for i in 1:5])
    f = R(1 // 2) * x[1]^7 * x[2]^2 + x[1] * x[2]^8 + x[3]^3
    @test QT.__chow_divided_difference(f, 1, 2) * (x[2] - x[1]) ==
      f - QT.__chow_swap(f, 1, 2)
    for e in 0:5
      f = prod(x[j] - x[i] for i in 1:e for j in (e + 1):5; init=R(1))
      monomials = QT.__chow_parabolic_monomials([5], [e])
      @test length(monomials) == binomial(5, e)
      @test all(
        QT.__chow_parabolic_alternant(f, prod(x .^ powers), [5], [e]) ==
        QT.__chow_divided_alternant(f * prod(x .^ powers), [5]) for powers in monomials
      )
    end
    @test length(QT.__chow_parabolic_monomials([5, 5, 2], [3, 3, 1])) == 200
  end

  @testset "Whole ideals, not just Hilbert functions" begin
    models = [
      (kronecker_quiver(4), [1, 1], [1, -1]),
      (kronecker_quiver(4), [1, 2], [2, -1]),
      (kronecker_quiver(3), [2, 3], [3, -2]),
      (subspace_quiver(5), [1, 1, 1, 1, 1, 2], [2, 2, 2, 2, 2, -5]),
      # A noncanonical, theta-coprime projective chamber; unsafe skips the
      # existing ample-stability API guard, not the Hilbert certificate.
      (subspace_quiver(5), [1, 1, 1, 1, 1, 2], [2, 2, 2, 2, 6, -7]),
      (Quiver([0 3 0; 0 0 0; 0 0 0]), [2, 3, 0], [3, -2, 0]),
    ]
    for (Q, d, theta) in models
      old = chow_ring(Q, d, theta; unsafe=true)[1][1]
      for adaptive in (false, true)
        result = QT.__chow_ring_parabolic(
          Q, d, theta, QT.extended_gcd(d)[2]; adaptive=adaptive
        )
        ring = result.presentation[1][1]
        @test same_ideal(old, ring)
        if adaptive
          @test result.certificate.graded_dimensions ==
            result.certificate.reineke_even_betti
        end
      end
      QT.empty_all_caches!()
    end
    Q = Quiver(zeros(Int, 1, 1))
    @test length(
      QT.gens(
        S.kbase(
          QT.quotient_ideal(chow_ring(Q, [1]; algorithm=:experimental_parabolic)[1][1])
        ),
      ),
    ) == 1
  end

  @testset "Normalization, cache, and published Chern controls" begin
    Q, d = kronecker_quiver(3), [2, 3]
    for chi in ([-1, 1], [2, -1])
      old = chow_ring(Q, d; chi=chi, unsafe=true)[1][1]
      new = chow_ring(Q, d; chi=chi, algorithm=:experimental_parabolic, unsafe=true)[1][1]
      @test same_ideal(old, new)
    end
    M = QuiverModuliSpace(Q, d)
    ring = chow_ring(M; algorithm=:experimental_parabolic)
    @test chow_ring(M) === ring
    @test same_ideal(chow_ring(Q, d; chi=[2, -1])[1][1],
      chow_ring(M; chi=[2, -1], algorithm=:experimental_parabolic))
    numbers = chern_numbers(M)
    @test numbers == Dict("c_1^6" => 41553, "c_2 c_1^4" => 21141,
      "c_2^2 c_1^2" => 10827, "c_2^3" => 5569, "c_3 c_1^3" => 6561,
      "c_3 c_2 c_1" => 3357, "c_3^2" => 1041, "c_4 c_1^2" => 1269,
      "c_4 c_2" => 643, "c_5 c_1" => 153, "c_6" => 13)
    M = QuiverModuliSpace(subspace_quiver(7), [ones(Int, 7); 2])
    chow_ring(M; algorithm=:experimental_parabolic, unsafe=true)
    @test chern_numbers(M; unsafe=true) == Dict("c_1^4" => 154, "c_2 c_1^2" => 112,
      "c_2^2" => 136, "c_3 c_1" => 56, "c_4" => 38)
  end

  @testset "Rank (5,5,2) sevenfolds" begin
    for adjacency in ([0 2 0; 0 0 1; 0 0 0], [0 2 1; 0 0 0; 0 0 0])
      Q, d = Quiver(adjacency), [5, 5, 2]
      @test !has_properly_semistable_representation(Q, d, canonical_stability(Q, d))
      result = QT.__chow_ring_parabolic(
        Q, d, canonical_stability(Q, d), QT.extended_gcd(d)[2]
      )
      @test result.certificate.graded_dimensions == [1, 2, 4, 4, 4, 4, 2, 1]
      @test result.certificate.schur_degree == 1
    end
  end

  @testset "Opt-in and errors" begin
    Q = kronecker_quiver(3)
    @test same_ideal(
      chow_ring(Q, [2, 3])[1][1], chow_ring(Q, [2, 3]; algorithm=:weyl)[1][1]
    )
    @test_throws ArgumentError chow_ring(Q, [2, 3]; algorithm=:unknown)
    @test_throws ArgumentError chow_ring(QuiverModuliSpace(Q, [2, 3]); algorithm=:unknown)
    @test_throws ArgumentError chow_ring(
      Q, [2, 3]; chi=[1, 1], algorithm=:experimental_parabolic
    )
    @test_throws ArgumentError chow_ring(Q, [2, 2]; algorithm=:experimental_parabolic)
    @test_throws ArgumentError chow_ring(Quiver([0 1; 1 0]), [1, 1], [1, -1];
      algorithm=:experimental_parabolic, unsafe=true)
  end
end
QT.empty_all_caches!()
end
