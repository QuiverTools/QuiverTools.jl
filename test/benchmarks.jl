#!usr/bin/env julia

using QuiverTools
using BenchmarkTools
using Memoization

setups = Dict{Function, Vector}()
results = Dict{Function, Any}()

"""bench_cached(f) runs the function f with the arguments in setups[f] and clears the cache of f afterwards.
It is used to benchmark functions that use memoization."""
function bench_cached(f)
	f(setups[f]...)
	Memoization.empty_cache!(f)
end

# TODO add benchmarks
setups[Euler_form] = [subspace_quiver(50), [[1 for i in 1:50]..., 50], [[-1 for i in 1:50]..., 1]]
setups[has_semistables] = [mKronecker_quiver(17), [7, 13], [13, -7]]
setups[all_HN_types] = [mKronecker_quiver(17), [7, 13], [13, -7]]



for f in keys(setups)
	results[f] = @benchmark bench_cached($f)
end

# TODO integrate results into the documentation