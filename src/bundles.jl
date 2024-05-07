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

using Combinatorics

part(n) = [x for x ∈ with_replacement_combinations(n:-1:0, n) if sum(x) == n]

partitions(n) = map(p -> filter(i -> i != 0, p), part(n))



# function wedge_of_tensor(p, q, n)
#     # this polynomial takes [λ^1(V),...,λ^{B.rank + r - 1}(V), λˆ1(W),...,λˆ{B.rank + r - 1}(W)] as formal variables
#     # and returns the polynomial form of \sum_{i,j} λ^i(V ⊠ W) λ^j(V ⊠ W) λ^{B.rank + r - 1 - i - j}(V ⊠ W). Products
#     # represent either tensor or box products, and sums represent direct sums.

#     λ(n)([],[])


function wedge_monomial(m::QQMPolyRingElem, n::Int)
    if length(terms(m)) != 1
        throw(ArgumentError("wedge_monomial only works for monomials"))
    end
    throw(ArgumentError("Not implemented yet"))
    # the problem here is that ∧(∧(things)) will appear, so a more general way to
    # compute exterior products is needed.
    # can I get away with computing all the polynomials P_{m,n} instead of just P_m?
    # can I do this with Schur functors?
end

function wedge(p::QQMPolyRingElem, n::Int)
    throw(ArgumentError("Not implemented yet"))
    if length(terms(p)) == 1
        return wedge_monomial(p, n)
    end
    return sum(wedge(first(p), i) * wedge(p - first(p), n - i) for i in 0:n)

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

    # compute the exterior powers of A
        # find a way to distribute ∧(sum(a_i))
        # apparently this can also be done with Schur functors.


    # compute the symmetric powers of B^\vee
        # apparently this can be done with Schur functors.

    # compute the terms of the Eagon-Northcott complex, as a polynomial not yet evaluated.

    # extract all the summands, only take them once.

    # only then evaluate the Teleman weights, if needed.
    # by keeping their polynomial description until the end like so,
    # we can also describe them in the Chow ring, to compute integrals.

end

