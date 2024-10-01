module BoundaryElementSolver

using LinearAlgebra, StaticArrays, Krylov, GaussQuadrature
using OhMyThreads: tmapreduce, @tasks
using Base.Threads: nthreads
using CairoMakie

# set multi-threading for BLAS
Krylov.BLAS.set_num_threads(nthreads())

export Model, Vertex, MeshedSurface, Triangle, PointCharge, PoissonSystem

export model_from_artifact, load_model
export sphere, sphere_surf
export Model2Surface, shift!

# high level functions
export solve, solve_ss, energy, surface_potential

export display_model, display_surface, display_system, display_potential


include("types.jl")
include("io.jl")
include("geometry.jl")
include("sphere.jl")

include("Poisson.jl")

include("utils/bisphere.jl")

include("utils/display.jl")

end