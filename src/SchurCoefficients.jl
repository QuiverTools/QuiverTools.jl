# Schur function coefficients used by polynomial semi-invariants.

const SIPartition = Tuple{Vararg{Int}}

_si_part(p::SIPartition, i::Int) = i <= length(p) ? p[i] : 0
_si_size(p::SIPartition) = sum(p; init=0)

function _si_trim(parts::Vector{Int})
  while !isempty(parts) && last(parts) == 0
    pop!(parts)
  end
  return Tuple(parts)
end

mutable struct _SchurContext
  rank::Int
  partitions::Dict{Int,Vector{SIPartition}}
  lr::Dict{Tuple{SIPartition,SIPartition,SIPartition},BigInt}
end

_SchurContext(rank::Int) = _SchurContext(
  rank, Dict{Int,Vector{SIPartition}}(),
  Dict{Tuple{SIPartition,SIPartition,SIPartition},BigInt}(),
)

function _si_partitions(context::_SchurContext, size::Int)
  return get!(context.partitions, size) do
    size == 0 && return SIPartition[()]
    return SIPartition[Tuple(p) for p in partitions(size) if length(p) <= context.rank]
  end
end

function _si_horizontal_strip(outer::SIPartition, inner::SIPartition)
  return all(
    _si_part(outer, i) >= _si_part(inner, i) >= _si_part(outer, i + 1)
    for i in 1:max(length(outer), length(inner))
  )
end

# Read the skew tableau right to left in each row, starting at the top.
function _si_count_tableaux(lambda::SIPartition, mu::SIPartition, nu::SIPartition)
  cells = [
    (row, col) for row in 1:length(nu)
    for col in nu[row]:-1:(_si_part(lambda, row) + 1)
  ]
  values = Dict{Tuple{Int,Int},Int}()
  used = zeros(Int, length(mu))

  function visit(position::Int)::BigInt
    position > length(cells) && return big(1)
    row, col = cells[position]
    right = get(values, (row, col + 1), length(mu))
    above = get(values, (row - 1, col), 0)
    count = big(0)
    for entry in (above + 1):right
      entry <= length(mu) || continue
      used[entry] < mu[entry] || continue
      used[entry] += 1
      if entry == 1 || used[entry - 1] >= used[entry]
        values[(row, col)] = entry
        count += visit(position + 1)
        delete!(values, (row, col))
      end
      used[entry] -= 1
    end
    return count
  end

  return visit(1)
end

function _si_lr_coefficient(
  context::_SchurContext, lambda::SIPartition, mu::SIPartition, nu::SIPartition
)
  return get!(context.lr, (lambda, mu, nu)) do
    _si_size(lambda) + _si_size(mu) == _si_size(nu) || return big(0)
    length(nu) <= context.rank || return big(0)
    all(_si_part(nu, i) >= _si_part(lambda, i) for i in 1:context.rank) ||
      return big(0)
    isempty(mu) && return lambda == nu ? big(1) : big(0)
    isempty(lambda) && return mu == nu ? big(1) : big(0)
    length(mu) == 1 && return _si_horizontal_strip(nu, lambda) ? big(1) : big(0)
    length(lambda) == 1 && return _si_horizontal_strip(nu, mu) ? big(1) : big(0)
    for (rectangle, other) in ((lambda, mu), (mu, lambda))
      if length(rectangle) == context.rank && all(==(rectangle[1]), rectangle)
        return if all(
          _si_part(nu, i) == rectangle[1] + _si_part(other, i)
          for i in 1:context.rank
        )
          big(1)
        else
          big(0)
        end
      end
    end
    return _si_count_tableaux(lambda, mu, nu)
  end
end

function _si_each_rectangle_partition(f, height::Int, width::Int)
  part = zeros(Int, height)
  function visit(row::Int, upper::Int)
    if row > height
      f(part)
      return nothing
    end
    for value in 0:upper
      part[row] = value
      visit(row + 1, value)
    end
  end
  visit(1, width)
end

# s_(w^h)^2 is multiplicity-free: the terms are (w+p, w-reverse(p))
# for partitions p inside the h-by-w rectangle.
function _si_rectangle_square(
  first::SIPartition, second::SIPartition, rank::Int
)
  first == second && !isempty(first) || return nothing
  all(==(first[1]), first) || return nothing
  height, width = length(first), first[1]
  product = Dict{SIPartition,BigInt}()
  _si_each_rectangle_partition(height, width) do part
    shape = _si_trim(
      vcat(
        [width + part[i] for i in 1:height],
        [width - part[i] for i in height:-1:1],
      ),
    )
    length(shape) <= rank && (product[shape] = big(1))
  end
  return product
end

function _si_schur_product(context::_SchurContext, factors)
  product = Dict{SIPartition,BigInt}(() => big(1))
  paired = nothing
  for i in eachindex(factors), j in (i + 1):length(factors)
    square = _si_rectangle_square(factors[i], factors[j], context.rank)
    if square !== nothing
      product, paired = square, (i, j)
      break
    end
  end
  for (index, factor) in enumerate(factors)
    paired !== nothing && index in paired && continue
    next = Dict{SIPartition,BigInt}()
    for (lambda, coefficient) in product
      for nu in _si_partitions(context, _si_size(lambda) + _si_size(factor))
        all(
          _si_part(nu, i) >= max(_si_part(lambda, i), _si_part(factor, i))
          for i in 1:context.rank
        ) || continue
        c = _si_lr_coefficient(context, lambda, factor, nu)
        c == 0 && continue
        next[nu] = get(next, nu, big(0)) + coefficient * c
      end
    end
    product = next
  end
  return product
end

function _si_schur_coefficient(context::_SchurContext, factors, target::SIPartition)
  isempty(factors) && return isempty(target) ? big(1) : big(0)
  length(factors) == 1 && return factors[1] == target ? big(1) : big(0)
  length(factors) == 2 && return _si_lr_coefficient(context, factors[1], factors[2], target)
  return get(_si_schur_product(context, factors), target, big(0))
end

function _si_complement(part::SIPartition, width::Int, rank::Int)
  length(part) <= rank && _si_part(part, 1) <= width || return nothing
  return _si_trim([width - _si_part(part, rank + 1 - i) for i in 1:rank])
end

# The coefficient of s_(width^rank) in a product of Schur functions.
function _si_rectangular_coefficient(context::_SchurContext, factors, width::Int)
  width < 0 && return big(0)
  target = width == 0 ? () : Tuple(fill(width, context.rank))
  isempty(factors) && return width == 0 ? big(1) : big(0)
  length(factors) == 1 && return factors[1] == target ? big(1) : big(0)
  length(factors) == 2 && return _si_complement(factors[1], width, context.rank) ==
         factors[2] ? big(1) : big(0)
  if length(factors) == 3
    complement = _si_complement(factors[3], width, context.rank)
    complement === nothing && return big(0)
    return _si_lr_coefficient(context, factors[1], factors[2], complement)
  end

  # Pair the two halves by complementary partitions. If one factor is a row,
  # remove it and use Pieri's horizontal-strip rule in the final pairing.
  row_index = findfirst(part -> length(part) == 1, factors)
  remaining = if row_index === nothing
    factors
  else
    [factors[i] for i in eachindex(factors) if i != row_index]
  end
  middle = length(remaining) ÷ 2
  left = _si_schur_product(context, remaining[1:middle])
  right = _si_schur_product(context, remaining[(middle + 1):end])
  result = big(0)
  for (lambda, a) in left
    complement = _si_complement(lambda, width, context.rank)
    complement === nothing && continue
    if row_index === nothing
      result += a * get(right, complement, big(0))
    else
      row_size = _si_size(factors[row_index])
      for (mu, b) in right
        _si_size(complement) - _si_size(mu) == row_size || continue
        _si_horizontal_strip(complement, mu) && (result += a * b)
      end
    end
  end
  return result
end
