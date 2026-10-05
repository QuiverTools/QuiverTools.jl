@testset "semi-invariant dimensions" begin
  Q = kronecker_quiver(3)
  @test semi_invariant_dimension(Q, [1, 1], [0, 0]) == 1
  @test semi_invariant_dimension(Q, [1, 1], [2, -2]) == 6
  @test semi_invariant_dimension(Q, [1, 1], [2, -2]; method=:lr) == 6
  @test semi_invariant_dimension(Q, [1, 1], [2, -2]; method=:reciprocity) == 6
  @test semi_invariant_dimension(Q, [1, 1], [-1, 1]) == 0
  @test semi_invariant_dimension(Q, [2, 1], [1, -2]) == 3
  @test semi_invariant_dimension(kronecker_quiver(1), [2, 2], [1, -1]) == 1
  @test semi_invariant_dimension(kronecker_quiver(1), [2, 1], [0, 0]) == 1
  @test semi_invariant_dimension(kronecker_quiver(2), [2, 1], [0, 0]) == 1
  # On the unframed two-arrow Kronecker quiver, the source coefficient
  # counts the six partitions inside a 2-by-2 rectangle.
  @test semi_invariant_dimension(kronecker_quiver(2), [2, 2], [2, -2]) == 6

  triangle = Quiver([0 1 1; 0 0 1; 0 0 0])
  @test semi_invariant_dimension(triangle, [1, 1, 1], [2, 0, -2]) == 3
  @test semi_invariant_dimension(triangle, [1, 0, 1], [1, 17, -1]) == 1
  path = Quiver([0 1 0; 0 0 1; 0 0 0])
  @test semi_invariant_dimension(path, [1, 2, 1], [1, 0, -1]) == 1
  @test semi_invariant_dimension(path, [1, 2, 1], [1, 0, -1];
    method=:reciprocity) == 1
  @test_throws DimensionMismatch semi_invariant_dimension(Q, [1], [1, -1])
  @test_throws ArgumentError semi_invariant_dimension(Q, [-1, 1], [1, -1])
  @test_throws ArgumentError semi_invariant_dimension(Quiver([0 1; 1 0]), [1, 1], [0, 0])

  # A nontrivial LR multiplicity, to check the tableau rule beyond Pieri.
  @test QuiverTools._si_lr_coefficient((2, 1), (2, 1), (3, 2, 1), 3,
    Dict{Any,BigInt}()) == 2
  square = QuiverTools._si_rectangle_square([(2, 2), (2, 2)], 4)
  @test length(square) == 6
  @test square[(4, 4)] == 1

  # The five rows of Table 4 in the extended-Dynkin manuscript.
  families = (
    (Quiver([0 2 0; 0 0 0; 0 2 0]), k -> [k, k, 1],
      (9, 75, 620, 5140)),
    (Quiver([0 2 0; 0 0 0; 1 1 0]), k -> [k, k, 1],
      (9, 75, 620, 5140)),
    (Quiver([0 1 1 0; 0 0 1 0; 0 0 0 0; 0 1 1 0]),
      k -> [k, k, k, 1], (8, 63, 504, 4090)),
    (Quiver([0 0 1 1 0; 0 0 1 1 0; 0 0 0 0 0;
        0 0 0 0 0; 0 0 1 1 0]),
      k -> [k, k, k, k, 1], (7, 52, 403, 3206)),
    (Quiver([0 0 0 0 1 0; 0 0 0 0 1 0; 0 0 0 0 1 0;
        0 0 0 0 1 0; 0 0 0 0 0 0; 0 0 0 0 1 0]),
      k -> [k, k, k, k, 2k, 1], (6, 42, 316, 2470)),
  )
  for (quiver, dimension, expected) in families, k in 1:4
    d = dimension(k)
    @test semi_invariant_dimension(quiver, d, canonical_stability(quiver, d)) ==
      expected[k]
  end
end
