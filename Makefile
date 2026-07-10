docs:
	julia --project=docs/ -e 'using Pkg; Pkg.update()'
	julia --project=docs/ docs/docs.jl
test:
	julia --project -e 'using Pkg; Pkg.activate("."); Pkg.test()'
format:
	julia .JuliaFormatter.jl
.PHONY: docs test format
