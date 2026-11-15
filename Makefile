THIS_BRANCH = $(shell git branch --show-current)

docs:
	# installs the last commit on this branch in the docs environment,
	#  replacing whatever version was there before
	cd docs && julia --project=. -e 'using Pkg; haskey(Pkg.project().dependencies, "QuiverTools") && Pkg.rm("QuiverTools"); Pkg.add(path="..", rev="$(THIS_BRANCH)"); Pkg.update()'
	julia --project=docs/ docs/docs.jl
test:
	julia --project -e 'using Pkg; Pkg.activate("."); Pkg.test()'
format:
	julia .JuliaFormatter.jl
.PHONY: docs test format
