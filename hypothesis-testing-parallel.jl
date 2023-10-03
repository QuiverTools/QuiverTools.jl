include("quivers.jl");

# Quivers = [GeneralizedKroneckerQuiver(i) for i in range(3,20)]
# d = [5,11]
# theta = [11,-5]

# #contains couples [dimension vector, stability parameter]
# parameters = [[[5,11],[11,-5]], [[7,9],[9,-7]],[[7,11],[11,-7]],[[5,12],[12,-5]],[[7,12],[12,-7]],[[5,13],[13,-5]],[[6,13],[13,-6]]]

# Quivers = [SubspaceQuiver(5)]
# d = [1,1,1,1,1,2]
# theta=[2,2,2,2,2,-5]

# parameters = [[d,theta]]

n = 10
Quivers = [ThreeVertexQuiver(i,j,k) for i in 1:n for j in 1:n for k in 1:n]
@info "number of quivers" length(Quivers)
dimension_list = [[i,j,k] for i in 1:n for j in 1:n for k in 1:n]
@info "number of dimensions" length(dimension_list)




Threads.@threads for Q in Quivers
        for d in dimension_list

# the following parallelization is worse for performance.
#  for Q in Quivers
#     Threads.@threads    for p in parameters
    theta = canonical_stability_parameter(Q,d)

    if has_semistable_representation(Q,d,theta) && is_coprime_for_stability_parameter(d,theta)
      
        HN = all_harder_narasimhan_types(Q, d, theta)
        AmplyStable = is_amply_stable(Q, d, theta)

        # technical condition: it implies AS, but is it equivalent?
        candidates_strong_AS = filter(e -> e != ZeroVector(number_of_vertices(Q)) && e != d && slope(e,theta) >= slope(d-e,theta),all_subdimension_vectors(d))
        strong_Ample_Stability = all(e -> euler_form(Q,e,d-e) <= -2, candidates_strong_AS)
        # technical condition: my condition on ext groups
#        MyHypothesis = all(hntype -> euler_form(Q,first(hntype),last(hntype)) < 0, filter(e -> length(e) > 1, HN))

        log = false
        if AmplyStable && !strong_Ample_Stability
            log = true
        end
        if log
            println(Q, " with d = ", d, " is not ok: ample stability is ", AmplyStable, " and strong ample stability is ", strong_Ample_Stability, ".")
        else
            #println(Q, " with d = ", d, " is ok.")
        end
    # else
    #     println("quiver: ", Q, " with dimension vector", d, " has no semistable representation")
    end

end 
end