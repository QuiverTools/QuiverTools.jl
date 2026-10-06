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
  @test QuiverTools._si_lr_coefficient(QuiverTools._SchurContext(3),
    (2, 1), (2, 1), (3, 2, 1)) == 2
  square = QuiverTools._si_rectangle_square((2, 2), (2, 2), 4)
  @test length(square) == 6
  @test square[(4, 4)] == 1

  # Multiple row factors at a sink.
  parallel_sources = Quiver([0 2 0; 0 0 0; 0 2 0])
  d = [2, 2, 1]
  @test semi_invariant_dimension(parallel_sources, d,
    canonical_stability(parallel_sources, d)) == 75

  # A mixed vertex with one incoming and two outgoing Schur factors.
  mixed_vertex = Quiver([0 2 0; 0 0 0; 1 1 0])
  @test semi_invariant_dimension(mixed_vertex, d,
    canonical_stability(mixed_vertex, d)) == 75

  # Four equal rectangular factors and one row at a sink.
  star = Quiver(
    [0 0 0 0 1 0; 0 0 0 0 1 0; 0 0 0 0 1 0;
      0 0 0 0 1 0; 0 0 0 0 0 0; 0 0 0 0 1 0],
  )
  d = [4, 4, 4, 4, 8, 1]
  @test semi_invariant_dimension(star, d,
    canonical_stability(star, d)) == 2470
end
