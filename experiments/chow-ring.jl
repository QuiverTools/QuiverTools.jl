using MegaQuiverTools, Nemo, IterTools, Oscar
# this wants to be an implementation of the Chow ring
# of an arbitrary quiver moduli, as presented by Hans
# in the 2013 Chow paper.


# TODO
# 2. Implement "minimal" forbidden subdimension vectors as in Hans's 2013 Chow paper
# 3. Compare outputs of rho with Sage code
# 4. Compare Itautological with Sage code

# data 
Q = GeneralizedKroneckerQuiver(3)
d = [2,3]
theta = CanonicalStabilityParameter(Q,d)
a = [-1,1]

# computations

R,x = polynomial_ring(Nemo.QQ, ["x$i$j" for i in 1:number_of_vertices(Q) for j in 1:d[i]])

# This is an ad-hoc way to index the varaibles.
function variable_index(d, i, j)
    return sum(d[1:i-1]) + j
end


function symmetric_polynomial(vars, degree::Int)
    return sum(prod(e) for e in IterTools.subsets(vars, degree))
end

# TODO change variable names to match Chow paper
A, e = polynomial_ring(Nemo.QQ, ["x$i$j" for i in 1:number_of_vertices(Q) for j in 1:d[i]])

targets = [[symmetric_polynomial([x[variable_index(d,i,j)] for j in 1:d[i]], k) for k in 1:d[i]] for i in 1:number_of_vertices(Q)]
targets = reduce(vcat, targets)

AinclusionR = hom(A,R,targets)

# forbidden = MegaQuiverTools.all_forbidden_subdimension_vectors(d,theta)

forbidden = [[1,1],[2,2]]

discriminant = prod(prod((x[variable_index(d,i,ell)] - x[variable_index(d,i,k)]) for k in 1:d[i] - 1 for ell in k+1:d[i]) for i in number_of_vertices(Q))


function symmetrisation(polynom, d, R, x)
    # this always returns zero. why?
    result = zero(R)
    for g in Iterators.product([SymmetricGroup(d[i]) for i in eachindex(d)]...)
        # build new index of variables
        newindex = [gi.d for gi in g]
        # adjust the index to account for the weird ordering of the variables
        if length(d) > 1
            for i in 2:length(d)
                newindex[i] = newindex[i] .+ sum(d[j] for j in 1:i-1)
            end
        end
        newindex = reduce(vcat, newindex)

        # build the new polynomial
        newvars = [x[newindex[i]] for i in eachindex(newindex)]
        result += prod(sign(gi) for gi in g) * polynom(newvars...)
    end
    return result
end

linear = sum(a[i] * e[variable_index(d,i,1)] for i in 1:number_of_vertices(Q))

Itautological = []

# nested loops bad but nested loops work
for dprime in forbidden
    for arrow in arrows(Q)
        for k in 1:dprime[arrow[1]]
            for ell in dprime[arrow[2]] + 1:d[arrow[2]]
                push!(Itautological, prod( (x[variable_index(d,arrow[2],ell)] - x[variable_index(d,arrow[1],k)]) for k in 1:dprime[arrow[1]] for ell in dprime[arrow[2]] + 1:d[arrow[2]]))
            end
        end
    end
end
@show Itautological
for taut in Itautological
    @show    symmetrisation(taut, d, R, sx)/discriminant
end
# Itautological = [prod( (x[variable_index(d,arrow[2],ell)] - x[variable_index(d,arrow[1],k)]) for k in 1:dprime[arrow[1]] for ell in dprime[arrow[2]] + 1:d[arrow[2]] for arrow in arrows(Q)) for dprime in forbidden]

kernel = [preimage(AinclusionR, symmetrisation(taut, d, R, x)/discriminant) + linear for taut in Itautological]