using MegaQuiverTools, IterTools, Nemo, Memoize

function bundles_on_stratum(Q::Quiver, hntype::Vector{Vector{Int}}, a::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    # weights
    kweights = map(di -> slope(di,theta, slope_denominator), hntype)
    kweights = kweights * lcm(denominator.(kweights))
    kweights = Int.(kweights)

    constant_term = sum(kweights[j]* (a' * hntype[j]) for j in eachindex(hntype))
    kweights = kweights .- constant_term

    U = Bundle[]
    for i in 1:nvertices(Q)
        # weight k_j appears d^j_i times
        new_weights = Vector{Int}(reduce(vcat, [[kweights[j] for k in 1:hntype[j][i] if hntype[j][i] > 0] for j in eachindex(hntype)]))
        new_rank = sum(hntype[j][i] for j in eachindex(hntype)) #this has to be d[i] by the way
        push!(U, Bundle(new_weights, new_rank))
    end
    return U
end

"""Returns a dictionary assigning to HN type (key) a Vector of Bundle objects.
Each Vector is an ordered list of U_i, and the weights are the Teleman weights of U_i on the HN type
given by the key."""
function all_universal_bundles(Q::Quiver, d::Vector{Int}, a::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    HN = filter(hntype -> hntype != [d], allHNtypes(Q, d, theta, slope_denominator))
    return Dict([hntype, bundles_on_stratum(Q, hntype, a, theta, slope_denominator)] for hntype in HN)
end

"""ONLY MEANT FOR THE 3-Kronecker QUIVER

Returns a list of terms in the generating family. If the keyword `left` is true, it returns
the left terms in the box product. Otherwise, it returns the right terms.
"""
function terms_from_symm(U::Vector{Bundle}, r::Int; left=true)

    if r == 1
        return [1]
    end

    # only relevant if r >= 2
    U1 = U[1]
    U2 = U[2]

    if left
        # k = 0 here
        kfirst = [⨂(MegaQuiverTools.dual(U2) for c in 1:r - 1)]
        middle = [⨂(MegaQuiverTools.dual(U1) for c in 1:k) ⊗ ⨂(MegaQuiverTools.dual(U2) for c in 1:r - 1 - k) for k in 1:r-2]
        # k = r-1 here
        klast = [⨂(MegaQuiverTools.dual(U1) for c in 1:r - 1)]
        out = vcat(kfirst, middle, klast)
    else
        kfirst = [⨂(U2 for c in 1:r - 1)]
        middle = [⨂(U1 for c in 1:k) ⊗ ⨂(U2 for c in 1:r - 1 - k) for k in 1:r-2]
        klast = [⨂(U1 for c in 1:r - 1)]
        out = vcat(kfirst, middle, klast)
    end
    return out
end



# this evaluates a Nemo monomial in a vector of Bundle objects

function evaluatemonomial(monomial, x, vars, left)
    ex = exponent_vector(monomial, 1)
    if left
        return ⨂(x[i] for i in filter(i -> ex[i] > 0, eachindex(x)) for c in 1:ex[i])
    else
        return ⨂(x[i] for i in filter(i -> ex[vars + i] > 0, eachindex(x)) for c in 1:ex[vars + i])
    end
end

"""ONLY MEANT FOR THE 3-Kronecker QUIVER

Returns a list part of the terms of the generating family coming from the wedge part of the r-th term of the Eagon-Northcott complex.
"""
function terms_from_wedge(U::Vector{Bundle}, A, B, r::Int; left=true)
    U1 = U[1]
    U2 = U[2]
    
    R,vars = polynomial_ring(QQ, vcat(["x$i" for i in 1:B.rank + r - 1], ["y$i" for i in 1:B.rank + r - 1]))
    
    @memoize Dict function polyλ(n)
        if n > 0 
            return 1//n * MegaQuiverTools.λ(n)(vars[1:n], vars[B.rank + r : B.rank + r + n])
        end
        return 1
    end
    
    # this polynomial takes [λ^1(V),...,λ^{B.rank + r - 1}(V), λˆ1(W),...,λˆ{B.rank + r - 1}(W)] as formal variables
    # and returns the polynomial form of \sum_{i,j} λ^i(V ⊠ W) λ^j(V ⊠ W) λ^{B.rank + r - 1 - i - j}(V ⊠ W). Products
    # represent either tensor or box products, and sums represent direct sums.

    P = sum(polyλ(i) * polyλ(j) * polyλ(B.rank + r - 1 - i - j) for i in 1:B.rank + r - 1 for j in 1:i if i + j <= B.rank + r - 1)
    # we don't actually need to multiply by 2 since we are interested in the monic monomials anyways.
    out = []
    
    for term in terms(P)
        coef = collect(coefficients(term))[1]
        mon = collect(monomials(term))[1]
        
        if left
            x = [wedge(U1, i) for i in 1:B.rank + r - 1]
        else
            x = [wedge(MegaQuiverTools.dual(U2), i) for i in 1:B.rank + r - 1]
        end
        
        # the coefficient determines the sign, i.e., if we take a dual or not.
        if coef > 0
            new = evaluatemonomial(mon, x, B.rank + r - 1, left)
        else
            new = dual(evaluatemonomial(mon, x, B.rank + r - 1, left))
        end
        
        if !(new in out)
            push!(out, new)
        end
        
        # the following adds terms with the correct multiplicity, should it be necessary.
        # for i in 1:abs(coef)
        #     push!(out, (Int(coef/abs(coef))) * evaluatemonomial(mon, x))
        # end

    end
    return out        
end

"""Returns the terms in the generating family coming from the r-th term of the Eagon-Northcott complex."""
function terms_from_r_term_in_EN_complex(U::Vector{Bundle}, A, B, r::Int; left=true)
    # terms from the wedge product
    wedge_terms = terms_from_wedge(U, A, B, r, left=left)
    # @info wedge_terms
    # terms from the symmetric product
    symm_terms = terms_from_symm(U, r, left=left)
    # @info symm_terms
    out = [x ⊗ y for x in wedge_terms for y in symm_terms]
    
    # broadcasting \otimes on Bundle objects does not seem to work directly. Too bad, for loop it is.
    if left
        factor = MegaQuiverTools.det(MegaQuiverTools.dual(U[1])) ⊗ MegaQuiverTools.det(MegaQuiverTools.dual(U[2]))
        out = map(elem -> elem ⊗ factor, out)
    else
        factor = MegaQuiverTools.det(U[1]) ⊗ MegaQuiverTools.det(U[2])
        out = map(elem -> elem ⊗ factor, out)
    end
    return out
end


"""Returns a list of terms in the generating family given the universal bundles U.
These depend on the HN stratum.

If the keyword `left` is true, it returns
the left terms in the box product. Otherwise, it returns the right terms.
"""
function full_family_on_stratum(Q::Quiver, U::Vector{Bundle}; left=true)
    A, B = MegaQuiverTools.KSmorphism(Q, U)
    out = []

    for r in 1:A.rank - B.rank + 1
        for new in terms_from_r_term_in_EN_complex(U, A, B, r, left=left)
            if !(new in out)
                push!(out, new)
            end
        end
    end
    # return reduce(vcat, terms_from_r_term_in_EN_complex(U, A, B, r, left=left) for r in 1:A.rank - B.rank + 1)
end
# for all r, the elements in the r term of the EN complex (on the given stratum)

function full_family(Q::Quiver, d::Vector{Int}, a::Vector{Int}, theta::Vector{Int}; slope_denominator::Function = sum, left=true)
    # throw(ArgumentError("Not implemented"))
    allbundles_allstrata = all_universal_bundles(Q, d, a, theta, slope_denominator)
    
    return Dict([hntype, full_family_on_stratum(Q, U, left=left)] for (hntype, U) in allbundles_allstrata)
    # for all r, the elements in the r term of the EN complex (on the given stratum)
end


"""Returns the "base change" function from the symmetric base to the power sum base for the ring of symmetric polynomials.
Denoting the symmetric polynomial base by ``e_i``, the power sum base by ``p_i`` and the base change function by ``\\nu_n``
such that ``p_n = \\nu_n(e_1,...,e_n)``, this function returns ``\\nu_n``."""
@memoize Dict function Newtonpolynomial(n)
    if n == 0
        throw(ArgumentError("Newtonpolynomial(0) is not defined"))
    elseif n == 1
        return x -> x[1]
    else
        function newPoly(x)
            return ((-1)^(n - 1) * n * x[n]) + reduce(+, (-1)^(i + n + 1) * x[n - i] * Newtonpolynomial(i)(x) for i in 1:n-1)
        end
    end
    return newPoly
end

"""Returns the function that takes two lists, [λ^i(V) for i in 1:n] and [λ^i(W) for i in 1:n],
and returns ``λ^n(V * W)`` as a function of ``λ^i(V)`` and ``λ^i(W)`` for ``i = 1,...,n``.
Note that for the notation to make sense, ``λ^1(E) = E`` for all ``E``.
"""
@memoize Dict function λ(n)
    if n == 0
        function λ0(x, y)
            return 1
        end
        return λ0
    elseif n == 1
        function λ1(x, y)
            return x[1] * y[1]
        end
        return λ1
    end
    function λn(x, y)
        out = Newtonpolynomial(n)(x) * Newtonpolynomial(n)(y)
        if n > 1
            out = out - reduce(+, (-1)^(n + i + 1) * λ(n - i)(x, y) * Newtonpolynomial(i)([λ(j)(x, y) for j in 1:i]) for i in 1:n-1)
        end
        return out
    end
    return λn
end

@memoize Dict function polyλ(n, vars)

    if n > 0 
        return 1//n * λ(n)(vars[1:n], vars[length(vars)/2 : length(vars)/2 + n])
    end
    return 1
end

function symm_of_tensor(factors, n)
    throw(ArgumentError("Not implemented"))
end


function wedge(p, n)
    throw(ArgumentError("Not implemented"))
    # pseudocode below


    # this should compute the wedge of a polynomial - that is, of a tensor product.
#    if p is a sum:
        # here p[i] is the i-th monomial    
#        return sum(wedge(p[1], i)*wedge(p[2:end], n - i) for i in 0:n)
#    else
        # here p is a monomial
#       if p is a product:
            # call the Schur functors on the first factor AND the rest of the factors.
#       else
            # return wedge of the monomial
#       end
#   end
end

function ENcomplex(Q::Quiver, U::Vector{Bundle})

    # this is going to be the lambda-ring
    # the n represents the n-th exterior power.
    variable_names = vcat(["u$i$n" for i in eachindex(U) for n in 1:U[i].rank], ["udual$i$n" for i in eachindex(U) for n in 1:U[i].rank])
    variable_names = map(x -> ["p$x", "q$x"], variable_names)
    variable_names = reduce(vcat, variable_names)

    # these are free variables, no relations between them have been established yet.
    R, vars = polynomial_ring(QQ, variable_names)

    # construct A and B^vee
    A = sum(Symbol("pu$i1")*Symbol("qudual$j1") for [i,j] in arrows(Q))
    Bdual = sum(Symbol("pudual$i1")*Symbol("qu$i1") for i in 1:nvertices(Q))
    # note that with this construction the rank can be found by evaluating the polynomial
    # with ui -> rank(U_i)
        
    # compute det(B^vee) with λn
    detBdual = wedge(Bdual, rank(Bdual))

    # compute the exterior powers of A
        # find a way to distribute ∧(sum(a_i))
        # apparently this can also be done with Schur functors.
    wedgesA = [wedge(A, i) for i in rank(B):rank(A)]


    # compute the symmetric powers of B^\vee
        # apparently this can be done with Schur functors.

    # compute the terms of the Eagon-Northcott complex, as a polynomial not yet evaluated.

    # extract all the summands, only take them once.

    # only then evaluate the Teleman weights, if needed.
    # by keeping their polynomial description until the end like so,
    # we can also describe them in the Chow ring, to compute integrals.

end


function main()

    Q = mKroneckerquiver(3)
    d = [2,3]
    a = [-1,1]
    theta = canonical_stability(Q, d)
    bounds = allTelemanbounds(Q, d, theta)

    @info "The computation of the full family is extremely slow."
    
    family = full_family(Q, d, a, theta, left=true)
    
    # test if objects have no self-extensions.
    exceptional_objects = true
    for hntype in keys(family)
        for elem in family[hntype]
            if elem != Bundle([],0) && maximum((dual(elem) ⊗ elem).weights) >= bounds[hntype]
                @info "Element $elem does not satisfy Teleman inequality on stratum $hntype."
                exceptional_objects = false
                break
            end            
        end
    end
    @info "\"All the objects are exceptional?\"" exceptional_objects

    # test if objects have no extensions between them.

    # test Teleman inequalities on members of the collection.
    # test chi vanishing of things in the collection.

end

main()