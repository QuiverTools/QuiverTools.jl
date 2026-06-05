using Test, QuiverTools, Documenter
using Pkg;
Pkg.activate(@__DIR__)

# Tests that have mathematical significance
# should be in the documentation doctests.

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

@testset "Universal Chern numbers: P^3 = M(kronecker_quiver(4), (1, 1))" begin
  M = QuiverModuliSpace(kronecker_quiver(4), [1, 1])
  chow_ring(M)
  cn = chern_numbers(M; universal=true)

  # the tangent-bundle Chern numbers are unchanged
  @test all(cn[k] == v for (k, v) in chern_numbers(M))

  # c_1(U_2) is linearly dependent on c_1(U_1) and is dropped; with the default
  # linearization c_1(U_1) = -h, so c_1(U_1)^3 = -1
  @test cn == Dict(
    "c_1^3" => 64,
    "c_2 c_1" => 24,
    "c_3" => 4,
    "c_1(U_1)^3" => -1,
  )
end;

@testset "Universal Chern numbers: 3-Kronecker, d = (2, 3)" begin
  M = QuiverModuliSpace(kronecker_quiver(3), [2, 3])
  chow_ring(M)
  cn = chern_numbers(M; universal=true)

  # the tangent-bundle Chern numbers are unchanged
  @test all(cn[k] == v for (k, v) in chern_numbers(M))

  # all intersection numbers in the universal Chern classes (the Chow ring
  # generators), with the default linearization chi = [-1, 1]. c_1(U_2) is
  # linearly dependent on c_1(U_1) and is dropped; these 14 match Meng's
  # Proposition 3.3.1.
  universal = Dict(k => v for (k, v) in cn if occursin("(U_", k))
  @test universal == Dict(
    "c_1(U_1) c_2(U_1) c_3(U_2)" => 2,
    "c_1(U_1) c_2(U_2) c_3(U_2)" => 3,
    "c_1(U_1)^2 c_2(U_1) c_2(U_2)" => 9,
    "c_1(U_1)^2 c_2(U_1)^2" => 6,
    "c_1(U_1)^2 c_2(U_2)^2" => 14,
    "c_1(U_1)^3 c_3(U_2)" => 5,
    "c_1(U_1)^4 c_2(U_1)" => 18,
    "c_1(U_1)^4 c_2(U_2)" => 27,
    "c_1(U_1)^6" => 57,
    "c_2(U_1) c_2(U_2)^2" => 5,
    "c_2(U_1)^2 c_2(U_2)" => 3,
    "c_2(U_1)^3" => 2,
    "c_2(U_2)^3" => 9,
    "c_3(U_2)^2" => 1,
  )
end;

@testset "HN types" begin
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

@testset "Constructors" begin
  # equivalences from the Sage docstrings; == compares adjacency matrices only
  @test thickened_subspace_quiver(2, 6) == three_vertex_quiver(0, 6, 6)
  @test generalized_subspace_quiver(2, [2, 3]) == three_vertex_quiver(0, 2, 3)
  @test generalized_subspace_quiver(3, [1, 1, 1]) == subspace_quiver(3)
  @test jordan_quiver() == loop_quiver(1)

  U = disjoint_union(kronecker_quiver(3), kronecker_quiver(4))
  @test n_vertices(U) == 4
  @test n_arrows(U) == 7
end;

@testset "string constructor" begin
  # chain semantics: a run of r hyphens is r arrows, chains are read left to right,
  # vertices are numbered in order of first appearance (parity with Sage's from_string)
  @test Quiver("a---b") == kronecker_quiver(3)
  @test Quiver("1--2-3") == Quiver([0 2 0; 0 0 1; 0 0 0])
  @test Quiver("a--b-3,a---3,3-a") == Quiver([0 2 3; 0 0 1; 1 0 0])
  @test Quiver("1--2,1---3,2----3") == Quiver([0 2 3; 0 0 4; 0 0 0])
end;
