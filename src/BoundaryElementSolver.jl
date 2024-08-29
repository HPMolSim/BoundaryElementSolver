module BoundaryElementSolver

using LinearAlgebra, StaticArrays, Krylov, GaussQuadrature
using GLMakie

export Model, Vertex, MeshedSurface, Triangle, PointCharge, PoissonSystem

export model_from_artifact, load_model
export sphere
export Model2Surface, shift!

# high level functions
export solve, energy, surface_potential

export display_model, display_surface, display_system, display_potential


include("types.jl")
include("io.jl")
include("geometry.jl")
include("sphere.jl")

include("Poisson.jl")
include("linear_solver.jl")

include("utils/bisphere.jl")

include("utils/display.jl")

end