format:
	julia .JuliaFormatter.jl

test:
	julia test.jl

.PHONY: format test
