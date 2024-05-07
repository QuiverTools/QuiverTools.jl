# This is a partial, naive implementation of tensor calculus.
# figure out a way to do this better.





# using IterTools

# import Base.==, Base.+, Base.-, Base.*, Base.show

# #this might be overengineered
# export ⊕, ⊗, ⊠, wedge, det, dual, symm, Bundle, EagonNorthcottcomplex, ⨂

# # for our purposes a Bundle is just a container for a set of weights and a rank.
# # TODO: implement a Zero bundle to streamline methods and cleanup the code.
# mutable struct Bundle
#     weights::Any # the terms here are supposed to be the weights in a weight space decomposition
#     rank::Int # this exists as a sanity check for implementation, as it should always be equal to length(weights)
    
#     # TODO: fix the Zero object behaviour and remove this. It bypasses the type check and that is prone to errors.
#     function Bundle(weights::Vector{Any}, rank::Int)
#         if length(weights) != rank
#             throw(ArgumentError("Length of weights must be equal to rank"))
#         end
#         return new(weights, rank)
#     end

#     function Bundle(weights::Vector{Int}, rank::Int)
#         if length(weights) != rank
#             throw(ArgumentError("Length of weights must be equal to rank"))
#         end
#         return new(weights, rank)
#     end
#     function Bundle(weights::Vector{Vector{Int}}, rank::Int)
#         if length(weights) != rank
#             throw(ArgumentError("Length of weights must be equal to rank"))
#         end
#         return new(weights, rank)
#     end
# end

# function ==(a::Bundle, b::Bundle)
#     return a.rank == b.rank && a.weights == b.weights
# end

# function show(io::IO, V::Bundle)
#     print(io, "Bundle of rank $(V.rank) with weights: $(V.weights)")
# end

# function dual(V::Bundle)
#     return Bundle(- V.weights, V.rank)
# end

# function ⊕(a::Bundle, b::Bundle)
#     return Bundle([a.weights; b.weights], a.rank + b.rank)
# end

# ###############################################
# # The following are useful for computations,  #
# # when 1 is treated as "the structure sheaf"  #
# # and 0 as "the zero bundle".                 #
# ###############################################

# function +(a::Bundle, b::Bundle)
#     return a ⊕ b
# end

# function *(a::Bundle, b::Bundle)
#     return a ⊗ b
# end

# function -(a::Bundle, b::Bundle)
#     return a ⊕ dual(b)
# end

# function *(a::Int, V::Bundle)
#     if a < 0
#         return dual((-a) * V)
#     end
#     return reduce(⊕, V for i in 1:a)
# end

# function *(V::Bundle, a::Int)
#     return a * V
# end

# function ⊗(b::Bundle, a::Int)
#     return a * b
# end
# function ⊗(a::Int, b::Bundle)
#     return a * b
# end

# function ⊗(a::Bundle, b::Bundle)
#     return Bundle([wa + wb for wa in a.weights for wb in b.weights], a.rank * b.rank)
# end

# function getweights(V::Bundle)
#     return V.weights
# end

# function ⨂(terms)
#     new = [sum(combination) for combination in IterTools.product(getweights.(terms)...)]
#     # the [:] thing is to flatten the array
#     return Bundle(new[:], prod(term.rank for term in terms))
# end

# function ⊠(a::Bundle, b::Bundle)
#     if a.weights isa Vector{Vector{Int}} || b.weights isa Vector{Vector{Int}}
#         throw(ArgumentError("Box products of box products are not implemented"))
#     end

#     return Bundle([[wa, wb] for wa in a.weights for wb in b.weights], a.rank * b.rank)
# end


# """Returns the k-th wedge power of the bundle V.

# Examples:

# ```julia-repl
# julia> U = Bundle([1, 2, 3], 3)
# Bundle of rank 3 with weights: [1, 2, 3]

# julia> wedge(U, 0)
# Bundle of rank 1 with weights: [0]
# julia> wedge(U, 1)
# Bundle of rank 3 with weights: [1, 2, 3]
# julia> wedge(U, 2)
# Bundle of rank 3 with weights: [3, 4, 5]
# julia> wedge(U, 3)
# Bundle of rank 1 with weights: [6]
# julia> wedge(U, 4)
# Bundle of rank 0 with weights: Int64[]
# ```
# """
# function wedge(V::Bundle, k::Int)
#     # TODO test the case k = 0 for box products properly
#     if k == 0 && V.weights isa Vector{Vector{Int}}
#         return Bundle([[0,0]],1)
#     end
#     if k > V.rank
#         # TODO simplify with ZeroBundle object
#         if V.weights isa Vector{Vector{Int}}
#             return Bundle(Vector{Int}[],0)
#         elseif V.weights isa Vector{Int}
#             return Bundle(Int[],0)
#         end
#     end
#     return Bundle([sum(e) for e in IterTools.subsets(V.weights,k)], binomial(V.rank,k))
# end


# """Returns the top exterior power of the bundle V."""
# function det(V::Bundle)
#     return wedge(V, V.rank)
# end


# function multisets_indices(n::Int, k::Int)
#     if k == 0
#         return [[]]
#     end
#     return [vcat([i], s) for i in 1:n for s in multisets_indices(n, k-1) if isempty(s) || i <= s[1]]
# end


# """Returns the k-th symmetric power of the bundle V.

# Examples:

# ```julia-repl
# julia> U = Bundle([1, 2, 3], 3);

# julia> symm(U,0)
# Bundle of rank 1 with weights: [0]
# julia> symm(U,1)
# Bundle of rank 3 with weights: [1, 2, 3]
# julia> symm(U,2)
# Bundle of rank 6 with weights: [2, 3, 4, 4, 5, 6]
# julia> symm(U,3)
# Bundle of rank 10 with weights: [3, 4, 5, 5, 6, 7, 6, 7, 8, 9]
# julia> symm(U,4)
# Bundle of rank 15 with weights: [4, 5, 6, 6, 7, 8, 7, 8, 9, 10, 8, 9, 10, 11, 12]
# ```
# """
# function symm(V::Bundle, k::Int)
#     if k == 0
#         if V.weights isa Vector{Vector{Int}}
#             return Bundle([[0, 0]], 1)
#         elseif V.weights isa Vector{Int}
#             return Bundle([0], 1)
#         end
#     end
#     return Bundle([sum(V.weights[i] for i in I) for I in multisets_indices(V.rank, k)], binomial(V.rank + k - 1, V.rank - 1))
# end

# # TODO test everything below this line

# # The methods below are the recursive method of finding λ^n(rs) as a function of λ^i(r), λ^i(s) for i = 0,…,n.
# # This is explained in
# # "Universal Polynomials in Lambda rings and the K-theory of the infinite loop space tmf", John R. Hopkinson,
# # Section 2 of https://dspace.mit.edu/bitstream/handle/1721.1/34544/71011847-MIT.pdf

# """Returns the "base change" function from the symmetric base to the power sum base for the ring of symmetric polynomials.
# Denoting the symmetric polynomial base by ``e_i``, the power sum base by ``p_i`` and the base change function by ``\\nu_n``
# such that ``p_n = \\nu_n(e_1,...,e_n)``, this function returns ``\\nu_n``."""
# @memoize Dict function Newtonpolynomial(n)
#     if n == 0
#         throw(ArgumentError("Newtonpolynomial(0) is not defined"))
#     elseif n == 1
#         return x -> x[1]
#     else
#         function newPoly(x)
#             return ((-1)^(n - 1) * n * x[n]) + reduce(+, (-1)^(i + n + 1) * x[n - i] * Newtonpolynomial(i)(x) for i in 1:n-1)
#         end
#     end
#     return newPoly
# end

# """Returns the function that takes two lists, [λ^i(V) for i in 1:n] and [λ^i(W) for i in 1:n],
# and returns ``λ^n(V * W)`` as a function of ``λ^i(V)`` and ``λ^i(W)`` for ``i = 1,...,n``.
# Note that for the notation to make sense, ``λ^1(E) = E`` for all ``E``.
# """
# @memoize Dict function λ(n)
#     if n == 0
#         function λ0(x, y)
#             return 1
#         end
#         return λ0
#     elseif n == 1
#         function λ1(x, y)
#             return x[1] * y[1]
#         end
#         return λ1
#     end
#     function λn(x, y)
#         out = Newtonpolynomial(n)(x) * Newtonpolynomial(n)(y)
#         if n > 1
#             out = out - reduce(+, (-1)^(n + i + 1) * λ(n - i)(x, y) * Newtonpolynomial(i)([λ(j)(x, y) for j in 1:i]) for i in 1:n-1)
#         end
#         return out
#     end
#     return λn
# end


# """Returns k times the exterior produt of ``A \\otimes B``.
# """
# function wedge_of_tensor(A::Bundle, B::Bundle, k)
#     out = Newtonpolynomial(k)([wedge(A, i) for i in 1:k]) ⊗ Newtonpolynomial(k)([wedge(B, i) for i in 1:k])
#     if k > 1
#         out = out - reduce(⊕, ((-1)^(k + i + 1) * wedge_of_tensor(A, B, k - i)) ⊗ Newtonpolynomial(i)([wedge_of_tensor(A, B, j) for j in 1:i]) for i in 1:k-1)
#     end
#     return out
# end
# """Returns k times the box produt of ``A \\boxtimes B``.
# """
# function wedge_of_box_product(A::Bundle, B::Bundle, k)
#     out = Newtonpolynomial(k)([wedge(A, i) for i in 1:k]) ⊠ Newtonpolynomial(k)([wedge(B, i) for i in 1:k])
#     if k > 1
#         out = out - reduce(⊕, ((-1)^(k + i + 1) * wedge_of_tensor(A, B, k - i)) ⊗ Newtonpolynomial(i)([wedge_of_tensor(A, B, j) for j in 1:i]) for i in 1:k-1)
#     end
#     return out
# end

# # Source and target of the Kodaira-Spencer morphism.
# # Input has to be a list of the summands of the universal representation.

# """
# Returns A and B from the standard 4-term exact sequence:

# `` 0 \\to
# \\mathit{Hom}(p^*\\mathcal{U},q^*\\mathcal{U}) \\to
# \\underbrace{\\bigoplus_{i \\in Q_0} \\mathcal{U}^\\vee_i \\boxtimes \\mathcal{U}_i}_{B^\\vee} \\to
# \\underbrace{\\bigoplus_{i \\to j \\in Q_1} \\mathcal{U}^\\vee_i \\boxtimes \\mathcal{U}_j}_{A^\\vee} \\to
# \\mathit{Ext}(p^*\\mathcal{U},q^*\\mathcal{U}) \\to
# 0.``
# """
# function KSmorphism(Q::Quiver, U::Vector{Bundle})
#     A = reduce(⊕, [U[a[1]] ⊠ dual(U[a[2]]) for a in arrows(Q)])
#     B = reduce(⊕, [U[i] ⊠ dual(U[i]) for i in 1:length(U)])
#     return [A,B]
# end


# """"Returns the r-th term of the Eagon-Northcott complex for the given source and target of the
# Kodaira--Spencer morphism.
# """
# function term_EagonNorthcottcomplex(A::Bundle, B::Bundle, r::Int)
#     # U is a list of the Ui, and r is the term we want to compute. In this convention,
#     # r=0 is the structure sheaf of the product, r=1 is the first term before that.
#     # r goes from 0 to rank(A) - rank(B) + 1.
    
#     if r == 0
#         return Bundle([[0,0]],1) #trivial bundle
#     elseif r > A.rank - B.rank + 1
#         throw(ArgumentError("r must be less than or equal to the rank of A"))
#     end
#     return wedge(A, B.rank + r - 1) ⊗ symm(dual(B),r - 1) ⊗ det(dual(B))
# end

# """Returns the terms of the Eagon-Northcott complex
# for the given quiver Q and the summands of the universal representation U.
# Excludes the ``0`` term, i.e., the trivial bundle.

# Examples:
# ```julia-repl
# julia> Q = mKroneckerquiver(2);
# julia> U = [Bundle([-1],1),Bundle([0],1)]
# 2-element Vector{Bundle}:
#  Bundle of rank 1, with weights: [-1]
#  Bundle of rank 1, with weights: [0]

# julia> EagonNorthcottcomplex(Q,U)
# 1-element Vector{Bundle}:
# Bundle of rank 1, with weights: [[-1, -1]]


# julia> Q = mKroneckerquiver(5);

# julia> EagonNorthcottcomplex(Q,U)
# 4-element Vector{Bundle}:
#  Bundle of rank 10, with weights: [[-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1]]
#  Bundle of rank 20, with weights: [[-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1]]
#  Bundle of rank 15, with weights: [[-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1]]
#  Bundle of rank 4, with weights: [[-1, -4], [-2, -3], [-3, -2], [-4, -1]]

# ```
# """
# function EagonNorthcottcomplex(Q::Quiver, U::Vector{Bundle})
#     AB = KSmorphism(Q, U)
#     return [term_EagonNorthcottcomplex(AB[1], AB[2], r) for r in 1:(AB[1].rank - AB[2].rank + 1)]
# end


