using MegaQuiverTools, JuMP, HiGHS

# Background: the Fano paper describes the ample cone in the case of ample stability; this thus can only be useful to find amples
# without ample stability.
# It is also a nice excuse to use linear programs.
"""Runs a linear program to find a linearization ``a`` of the universal bundles whose higher cohomology vanishes."""
function linearization_vanishing_universal_bundles(Q,d,theta)
    
    if gcd(d...) != 1
        throw(ArgumentError("The dimension vector must be coprime, otherwise universal bundles don't exist"))
    elseif !MegaQuiverTools.HasSemistableRepresentation(Q,d,theta)
        throw(ArgumentError("The moduli space is empty."))
    end
    
    HN = filter(dstar -> length(dstar) > 1, AllHarderNarasimhanTypes(Q,d,theta))
    
    eta = MegaQuiverTools.all_weight_bounds(Q,d,theta)
    
    K1 = []
    Arows = []
    
    for dstar in HN
        kappas = map(di -> Slope(di,theta),dstar)
        kappas = kappas * lcm(denominator.(kappas))
        push!(K1,kappas[1])
        
        D = hcat(dstar...)
        push!(Arows, D*kappas)
    end
    
    A = hcat(Arows...)'
    
    m = Model(HiGHS.Optimizer)
    @variable(m, x[1:number_of_vertices(Q)],integer = true)
    @constraint(m, A*x .<= eta .+ K1 .- 1)
    @constraint(m, x' * d == 1)
    @objective(m, Min, 0)
    set_silent(m)
    optimize!(m)
    
    if JuMP.is_solved_and_feasible(m)
        return value.(x)
    else
        return nothing
    end
    
end


function main()

    # special cases
    qlist = [ThreeVertexQuiver(1,6,7)]
    dlist = [[4,6,7]]

    @info qlist[1]
    @info dlist[1]
    @info CanonicalStabilityParameter(qlist[1],dlist[1])
    @info MegaQuiverTools.is_schur_root(qlist[1],dlist[1])
    @info IsAmplyStable(qlist[1],dlist[1],CanonicalStabilityParameter(qlist[1],dlist[1]))
    @info MegaQuiverTools.is_coprime_for_stability_parameter(dlist[1],CanonicalStabilityParameter(qlist[1],dlist[1]))

    # Three vertex quivers
    # qlist = [ThreeVertexQuiver(i,j,k) for i in 1:10 for j in 1:10 for k in 1:10]
    # dlist = [[i,j,k] for i in 1:10 for j in 1:10 for k in 1:10]

    # Generalized Kronecker quivers
    # qlist = [GeneralizedKroneckerQuiver(n) for n in 1:10]
    # dlist = [[n,m] for n in 1:10 for m in 1:10]

    for Q in qlist
        for d in dlist
            theta = CanonicalStabilityParameter(Q,d)
            # ensure that universal bundles exist and the moduli space is non-empty
            if gcd(d) == 1 && HasSemistableRepresentation(Q,d,theta)
                x = linearization_vanishing_universal_bundles(Q,d,theta)
                @info x
            end
        end
    end

end


main()
