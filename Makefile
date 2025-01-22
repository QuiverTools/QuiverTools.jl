docs:
	julia --project=docs/ docs/docs.jl
test:
	julia --project -e 'using Pkg; Pkg.activate("."); Pkg.test()'
format:
	julia .JuliaFormatter.jl
.PHONY: docs test format