#!usr/bin/env julia

using QuiverTools
using BenchmarkTools
using Memoization
using Serialization

BASE_PATH = dirname(dirname(pathof(QuiverTools)))

SUITE = BenchmarkGroup()
SUITE[:benchmark] = BenchmarkGroup()

"""bench_cached(f) runs the function f with the arguments in setups[f] and clears the cache of f afterwards.
It is used to benchmark functions that use memoization."""
function bench_cached(f, args...)
	f(args...)
	Memoization.empty_all_caches!()
end

@info "Setting up functions and data for benchmarks..."

batch_0_data = [
	[mKronecker_quiver(5)],
	# [subspace_quiver(20)],
	# [three_vertex_quiver(1, 6, 7)]
];

quiver_data_light = [
	(mKronecker_quiver(5), [2, 3], [3, -2]),
	# (subspace_quiver(5), [[1 for i in 1:5]..., 5], [[-1 for i in 1:5]..., 1]),
	# (three_vertex_quiver(1, 2, 1), [2, 1, 3], [1, 1, -1]),
];

# quiver_data_heavy = [
# 	(mKronecker_quiver(17), [7, 13], [13, -7]),
# 	(subspace_quiver(10), [[1 for i in 1:10]..., 10], [[-1 for i in 1:10]..., 1]),
# 	(three_vertex_quiver(1, 6, 7), [4, 1, 4], [25, 24, -31]),
# ];

quiver_data_heavy = [
	(mKronecker_quiver(3), [2, 3], [3, -2]),
];

functions_batch_0 = [
	is_acyclic,
	is_connected,
];

functions_batch_1 = [
	all_HN_types,
	all_Teleman_bounds,
	Euler_form,
	has_semistables,
];

functions_batch_2 = [
	QuiverTools._Hodge_polynomial_fast,
	# Hodge_polynomial,
];

functions_batch_3 = [
	# generic_ext,
	# generic_hom,
	is_Schur_root,
	canonical_decomposition,
	in_fundamental_domain,
];

setups = Dict()


for f in functions_batch_0
	setups[f] = batch_0_data
end

for f in functions_batch_1
	setups[f] = vcat(quiver_data_heavy, quiver_data_light)
end

for f in functions_batch_2
	setups[f] = quiver_data_light
end

for f in functions_batch_3
	setups[f] = map(t -> (t[1], t[2]), vcat(quiver_data_light, quiver_data_heavy))
end

# TODO add benchmarks

functions = vcat(functions_batch_0, functions_batch_1, functions_batch_2, functions_batch_3);

SUITES = BenchmarkGroup()

@info "Setting up benchmarks..."
for f in functions
	SUITES[f] = BenchmarkGroup()
	for datum in setups[f]
		SUITES[f][datum] = @benchmarkable bench_cached($f, $datum...)
	end
end
@info "Tuning benchmarks..."
tune!(SUITES)

@info "Running benchmarks..."
res = run(SUITES)

@info "Writing results to file..."

INITIAL_MESSAGE = "These are the benchmarks for the functions in QuiverTools.\n\nThe benchmarks were generated by running\n\n```shell\njulia --project benchmark/benchmarks.jl\n```\n\nfrom the root of the repository.\nThis was done on a 2021 MacBook Pro with a M1 Pro chip and 16 GB of RAM.\n\n";

open(BASE_PATH * "/docs/src/benchmarks.md", "w") do file
	write(file, "# Benchmarks\n\n")
	write(file, INITIAL_MESSAGE)
	for f in functions
		write(file, "## Benchmarking `$f`\n\n")
		for datum in setups[f]
			write(file, "### Datum: $datum\n\n```julia\n")
			show(file, "text/plain", res[f][datum])
			write(file, "\n```\n")
		end
	end
end

for f in functions
	for datum in setups[f]
		open(BASE_PATH * "/benchmark/raw/$(f)_$(datum).jls", "w") do file
			serialize(file, res[f][datum])
		end
	end
end
