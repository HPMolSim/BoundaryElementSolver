JL = julia --project

default: init update

init:
	$(JL) -e 'using Pkg; Pkg.develop([PackageSpec(; path="..")]); Pkg.instantiate(); Pkg.precompile();'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile();'