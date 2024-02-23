struct MeshedSurface{T}
    vertices::Vector{SVector{3, T}}
    faces::Vector{SVector{3, Int}}
    points::Vector{SVector{3, T}}
    normals::Vector{SVector{3, T}}
    areas::Vector{T}
end

struct Charge{T}
    position::SVector{3, T}
    charge::T
    ϵ::T
end

struct PoissonSystem{T}
    N_s::Int # number of surfaces
    N_c::Int # number of charges
    N_m::Vector{Int} # number of points on each surface

    surfaces::Vector{MeshedSurface{T}}
    centers::Vector{SVector{3, T}}
    
    ϵ_medium::T
    ϵ_particle::Vector{T}

    charges::Vector{Charge{T}}
end

struct PoissonSolver{T}
    Lhs::AbstractArray{T, 2}
    Rhs::AbstractArray{T, 1}
    Preconditioner::AbstractArray{T, 2}
    max_error::T
end