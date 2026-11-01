# Run with `julia --project examples/quasimap_quantum_cohomology.jl`.
# All computations are exact, including the characteristic polynomials.
using QuiverTools

const N = QuiverTools.Singular.Nemo
const T, t = N.polynomial_ring(N.QQ, "t")

for (n, d) in ((3, [1, 1]), (4, [1, 2]), (3, [2, 3]), (4, [2, 3]))
  M = QuiverModuliSpace(kronecker_quiver(n), d)
  A = quasimap_quantum_cohomology(M)
  x = quantum_tautological_classes(A)
  # In this positive regime the degree-one quantum classes are classical.
  c1 = -sum(canonical_stability(M.Q, d)[i] * x[i][1] for i in eachindex(d))
  println("\n", n, "-Kronecker moduli, dimension vector ", d)
  println("  Gromov--Witten comparison: ", gw_comparison(A))
  matrix, basis, _ = quantum_multiplication_matrix(A, c1, [1])
  println("  Rank at q=1: ", length(basis))
  println("  Characteristic polynomial of c1*: ", N.charpoly(T, matrix))
  classical, _, _ = quantum_multiplication_matrix(A, c1, [0])
  println("  Characteristic polynomial at q=0: ", N.charpoly(T, classical))
end

# Multiple Novikov parameters: P^1 x P^1.
A = quasimap_quantum_cohomology(
  QuiverModuliSpace(three_vertex_quiver(2, 2, 0), [1, 1, 1]); chi=[1, 0, 0]
)
x = quantum_tautological_classes(A)
matrix, basis, _ = quantum_multiplication_matrix(A, x[2][1] + x[3][1], [1, 1])
println("\nP^1 x P^1, hyperplane sum at q1=q2=1: ", N.charpoly(T, matrix))
