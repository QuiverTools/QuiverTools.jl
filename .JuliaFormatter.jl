using JuliaFormatter

# every Julia file that is tracked, or untracked but not ignored by .gitignore
files = readlines(`git ls-files --cached --others --exclude-standard -- '*.jl'`)

clean = map(format, files)
all(clean) && exit(0)

# print non-formatted files
for x in files[.!(clean .== true)]
  println("Not formatted: ", x)
end
exit(1)
