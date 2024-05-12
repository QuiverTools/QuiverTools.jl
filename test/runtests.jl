using Test, QuiverTools

@testset "basic-methods" begin
    K = mKronecker_quiver(4)
    
    @test QuiverTools.underlying_graph(K) == [0 4;4 0]
    @test nvertices(K) == 2
    @test narrows(K) == 4
    @test is_acyclic(K) == true
    @test is_connected(K) == true
    @test is_source(K, 1) == true && is_source(K, 2) == false
    @test is_sink(K, 2) == true && is_sink(K, 1) == false
    @test QuiverTools.Euler_matrix(K) == [1 -4; 0 1]
    @test Euler_form(K, [1, 1], [1, 1]) == -2
end;

@testset "constructors" begin
    @test mKronecker_quiver(3).adjacency == [0 3; 0 0]
    @test three_vertex_quiver(3,4,5).adjacency == [0 3 4; 0 0 5; 0 0 0]
    @test loop_quiver(3).adjacency == Matrix{Int64}(reshape([3],1,1))
    @test subspace_quiver(2).adjacency == [0 0 1; 0 0 1;0 0 0]
end;

@testset "all_slope_decreasing_sequences()" begin
    Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];
    @test QuiverTools.all_slope_decreasing_sequences(Q, d, theta) == [  [[2, 3]],
                                                                        [[1, 1], [1, 2]],
                                                                        [[2, 2], [0, 1]],
                                                                        [[2, 1], [0, 2]],
                                                                        [[1, 0], [1, 3]],
                                                                        [[1, 0], [1, 2], [0, 1]],
                                                                        [[1, 0], [1, 1], [0, 2]],
                                                                        [[2, 0], [0, 3]]]
end;

@testset "all_generic_subdimension_vectors()" begin
    Q = mKronecker_quiver(3);
    @test QuiverTools.all_generic_subdimension_vectors(Q, [2,3]) == [   [0, 0],
                                                                        [0, 1],
                                                                        [0, 2],
                                                                        [1, 2],
                                                                        [0, 3],
                                                                        [1, 3],
                                                                        [2, 3]]
    @test QuiverTools.all_generic_subdimension_vectors(Q, [3,0]) == [   [0, 0],
                                                                        [1, 0],
                                                                        [2, 0],
                                                                        [3, 0]]
end;

@testset "has_semistables()" begin
    Q = mKronecker_quiver(17)
    @test has_semistables(Q,[1,13], [13,-1]) == true
    @test has_semistables(Q,[1,13], [-13,1]) == false

    K = mKronecker_quiver(4)
    @test has_semistables(K,[1,5], [5,-1]) == false
    @test has_semistables(K,[1,5], [-5,1]) == false

    A2 = mKronecker_quiver(1); d = [1,1]; theta = [1,-1]
    @test has_semistables(A2, d, theta) == true

    d = [2,2]
    @test has_semistables(A2, d, theta) == true

    d = [1,2]
    @test has_semistables(A2, d, theta) == false

    d = [0,0]
    @test has_semistables(A2, d, theta) == true


    K3 = mKronecker_quiver(3); theta = [3,-2]; d = [2,3];
    @test has_semistables(K3, d, theta) == true

    d = [1,4]
    @test has_semistables(K3, d, theta) == false

    @test has_semistables(K3, [3,0], [0,-1]) == true
    @test has_semistables(K3, [0,3], [1,0]) == true
end;

@testset "has_stables()" begin
    Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];
    @test has_stables(Q, d, theta) == true
    
    @test has_stables(Q, [1,0], [0,-1]) == true
    @test has_stables(Q, [0,1], [1,0]) == true
    @test has_stables(Q, [3,0], [0,-1]) == false
    @test has_stables(Q, [0,3], [1,0]) == false
end

@testset "proper-semistability" begin
    Q = mKronecker_quiver(2); d = [2,2]; theta = [1,-1];

    @test has_semistables(Q, d, theta) == true
    @test has_stables(Q, d, theta) == false
    Q = mKronecker_quiver(3);
    @test has_stables(Q, [1,0], [0,-1]) == true
    @test has_stables(Q, [0,1], [1,0]) == true
    @test has_stables(Q, [3,0], [0,-1]) == false
    @test has_stables(Q, [0,3], [1,0]) == false
    @test has_semistables(Q, [3,0], [0,-1]) == true
    @test has_semistables(Q, [0,3], [1,0]) == true
end;
@testset "all_HN_types()" begin
    Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

    @test all_HN_types(Q, d, theta) == [[[2, 3]],
                                        [[1, 1], [1, 2]],
                                        [[2, 2], [0, 1]],
                                        [[2, 1], [0, 2]],
                                        [[1, 0], [1, 3]],
                                        [[1, 0], [1, 2], [0, 1]],
                                        [[1, 0], [1, 1], [0, 2]],
                                        [[2, 0], [0, 3]]]

    theta = [-3,2]
    @test all_HN_types(Q, d, theta) == [[[0,3],[2,0]]]

    Q = three_vertex_quiver(3,4,5); d = [3,5,7]; theta = [43,26,-37];
    
    # 3vertexquiver-3-5-7-canonical.txt
    expected = "" #has to be initialised outside of the open file
    open("3vertexquiver-3-5-7-canonical.txt","r") do file
        expected = readline(file)
    end

    @test string(all_HN_types(Q,d,theta)) == expected
end;

@testset "is_amply_stable()" begin
    Q = mKronecker_quiver(3); d = [2,3];

    @test is_amply_stable(Q, d, [3,-2]) == true
    @test is_amply_stable(Q, d, [-3,2]) == false
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