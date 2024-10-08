using Test, QuiverTools, Documenter

@info "Almost all the tests are in the documentation."

DocMeta.setdocmeta!(QuiverTools, :DocTestSetup, :(using QuiverTools))
doctest(QuiverTools, manual=false, testset="Doctests")

@testset "strict sst" begin
    # proper-semistability
    Q = mKronecker_quiver(2)
    d = [2, 2]
    theta = [1, -1]

    @test has_semistables(Q, d, theta) == true
    @test has_stables(Q, d, theta) == false
    Q = mKronecker_quiver(3)
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

    @test string(all_HN_types(Q, d, theta)) == expected
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
