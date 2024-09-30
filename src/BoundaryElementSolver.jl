module BoundaryElementSolver

using LinearAlgebra, StaticArrays, Krylov, GaussQuadrature
using CairoMakie

export Model, Vertex, MeshedSurface, Triangle, PointCharge, PoissonSystem

export model_from_artifact, load_model
export sphere, sphere_surf
export Model2Surface, shift!

# high level functions
export solve, energy, surface_potential

export display_model, display_surface, display_system, display_potential


include("types.jl")
include("io.jl")
include("geometry.jl")
include("sphere.jl")

include("Poisson.jl")

include("utils/bisphere.jl")

include("utils/display.jl")

end