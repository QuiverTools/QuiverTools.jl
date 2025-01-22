docs:
	julia --project=docs/ docs/docs.jl
test:
	julia --project test/runtests.jl
format:
	julia .JuliaFormatter.jl
.PHONY: docs test format
