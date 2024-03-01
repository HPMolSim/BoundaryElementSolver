# vertices, normals and faces will be loaded from files.
struct Model{T}
    verts::Dict{Int64, Tuple{SVector{3, T}, SVector{3, T}}}
    faces::Vector{SVector{3, Int}}
end

struct Triangle{T}
    r::SVector{3, T}
    n::SVector{3, T}
    a::T
end

mutable struct MeshedObject{T}
    # info of surfaces
    tris::Vector{Triangle{T}}
    center::SVector{3, T}
    angle::SVector{3, T}
    ϵ::T
end

mutable struct PointCharge{T}
    position::SVector{3, T}
    charge::T
    ϵ::T
end

mutable struct PoissonSystem{T}
    ϵ_medium::T
    N_s::Int # number of surfaces
    N_c::Int # number of charges
    N_m::Vector{Int} # number of points on each surface

    surfaces::Vector{MeshedObject{T}}
    charges::Vector{PointCharge{T}}
end

mutable struct PoissonSolver{T}
    Lhs::AbstractArray{T, 2}
    Rhs::AbstractArray{T, 1}
    Preconditioner::AbstractArray{T, 2}
    max_error::T
end