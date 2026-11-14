# Exact ideal comparisons outside the timing runs. The censored large cases
# are deliberately absent: no reference presentation was obtained for them.
using QuiverTools
const QT = QuiverTools
const S = QT.Singular
cases = [
  ("P3", kronecker_quiver(4), [1, 1]),
  ("Gr24", kronecker_quiver(4), [1, 2]),
  ("Gr36", kronecker_quiver(6), [1, 3]),
  ("K3_23", kronecker_quiver(3), [2, 3]),
  ("subspace7", subspace_quiver(7), [ones(Int, 7); 2]),
  ("K3_34", kronecker_quiver(3), [3, 4]),
]
for (name, Q, d) in cases
  first = chow_ring(Q, d; algorithm=:weyl, unsafe=true)[1][1]
  second = chow_ring(Q, d; algorithm=:experimental_parabolic, unsafe=true)[1][1]
  I, J = QT.quotient_ideal(first), QT.quotient_ideal(second)
  A, B = QT.base_ring(I), QT.base_ring(J)
  to_second = S.AlgebraHomomorphism(A, B, QT.gens(B))
  to_first = S.AlgebraHomomorphism(B, A, QT.gens(A))
  @assert all(iszero(S.reduce(to_second(f), J)) for f in QT.gens(I))
  @assert all(iszero(S.reduce(to_first(f), I)) for f in QT.gens(J))
  println(name, ": whole-ideal equality passed")
  QT.empty_all_caches!()
end
