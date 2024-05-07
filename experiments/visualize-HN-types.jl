using QuiverTools, Plots

"""Converts a list of Harder-Narasimhan types to a list of tuples of the form [zero, sum(hntype[1]), sum(hntype[1:2]), ...] in tuples."""
function reformat(hntype)
    out = [Tuple((0 for i in 1:length(hntype[1])))]
    for i in 1:length(hntype)
        push!(out, Tuple(sum(hntype[j] for j in 1:i)))
    end
    return out
end



myplot = plot()
for Q in [ThreeVertexQuiver(1,6,7)]
    for d in [[4,6,7]]
        theta = CanonicalStabilityParameter(Q,d)
        hntypes = AllHarderNarasimhanTypes(Q,d,theta)
        for hntype in hntypes
            hntype = reformat(hntype)
            plot!(myplot, hntype, label = "HN type")
        end
    end
end