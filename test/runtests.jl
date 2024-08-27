using BoundaryElementSolver
using LinearAlgebra, StaticArrays
using Test

@testset "BoundaryElementSolver.jl" begin
    include("types.jl")
    include("sphere.jl")
    include("geometry.jl")
    include("Poisson.jl")
end
