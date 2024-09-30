# vertices, normals and faces will be loaded from files.
struct Vertex{T}
    p::SVector{3, T} # position
    n::SVector{3, T} # normal
end

Base.show(io::IO, v::Vertex{T}) where{T} = print(io, "Vertex{$T}, p: $(v.p), n: $(v.n)")
Base.:(==)(v1::Vertex{T}, v2::Vertex{T}) where{T} = v1.p == v2.p && v1.n == v2.n

Vertex(p::NTuple{3, T}, n::NTuple{3, T}) where{T} = Vertex(SVector{3, T}(p), SVector{3, T}(n ./ norm(n)))
Base.copy(v::Vertex{T}) where{T} = Vertex{T}(copy(v.p), copy(v.n))
mid(a::Vertex{T}, b::Vertex{T}) where{T} = Vertex{T}(SVector{3, T}((a.p + b.p)/2), SVector{3, T}((a.n + b.n)/norm(a.n + b.n)))
mid(a::Vertex{T}, b::Vertex{T}, c::Vertex{T}) where{T} = Vertex{T}(SVector{3, T}((a.p + b.p + c.p)/3), SVector{3, T}((a.n + b.n + c.n)/norm((a.n + b.n + c.n))))

struct Model{T}
    verts::Vector{Vertex{T}}
    faces::Vector{SVector{3, Int}}
end

Base.show(io::IO, v::Model{T}) where{T} = print(io, "Model{$T}, nv: $(length(v.verts)), nf: $(length(v.faces))")
Base.copy(m::Model{T}) where{T} = Model{T}(copy(m.verts), copy(m.faces))


mutable struct Triangle{T}
    r::SVector{3, T}
    n::SVector{3, T}
    a::T

    function Triangle(v1::Vertex{T}, v2::Vertex{T}, v3::Vertex{T}) where{T}
        r = (v1.p + v2.p + v3.p) ./ 3
        n = (v1.n + v2.n + v3.n) ./ 3
        n = n ./ norm(n)
        a = tri_area(v1.p, v2.p, v3.p)
        return new{T}(r, n, a)
    end
end

Base.show(io::IO, v::Triangle{T}) where{T} = print(io, "Triangle{$T}, r: $(v.r), n: $(v.n), a: $(v.a)")
Base.copy(t::Triangle{T}) where{T} = Triangle{T}(copy(t.r), copy(t.n), t.a)

mutable struct MeshedSurface{T}
    tris::Vector{Triangle{T}} # info of surfaces
end

Base.show(io::IO, v::MeshedSurface{T}) where{T} = print(io, "MeshedSurface{$T}, n: $(length(v.tris))")
nt(m::MeshedSurface) = length(m.tris)
Base.copy(m::MeshedSurface{T}) where{T} = MeshedSurface{T}(copy(m.tris))

struct PointCharge{T}
    r::SVector{3, T}
    q::T
    function PointCharge(x::T, y::T, z::T, q::T) where{T}
        new{T}(SVector{3, T}(x, y, z), q)
    end
    function PointCharge(r::SVector{3, T}, q::T) where{T}
        new{T}(r, q)
    end
end

Base.show(io::IO, v::PointCharge{T}) where{T} = print(io, "PointCharge{$T}, position: $(v.r), charge: $(v.q)")
Base.copy(p::PointCharge{T}) where{T} = PointCharge(copy(p.r), p.q)

# T for the type of the positions
# TE for eps, can be complex numbers
struct PoissonSystem{T, TE}
    ϵ_medium::TE

    ϵ_surfaces::Vector{TE}
    surfaces::Vector{MeshedSurface{T}}

    ϵ_charges::Vector{TE}
    charges::Vector{PointCharge{T}}

    function PoissonSystem(ϵ_medium::TE, ϵ_surfaces::Vector{TE}, surfaces::Vector{MeshedSurface{T}}, ϵ_charges::Vector{TE}, charges::Vector{PointCharge{T}}) where{T, TE}
        if length(ϵ_surfaces) != length(surfaces)
            throw(ArgumentError("length(ϵ_surfaces) != length(surfaces)"))
        end
        if length(ϵ_charges) != length(charges)
            throw(ArgumentError("length(ϵ_charges) != length(charges)"))
        end
        new{T, TE}(ϵ_medium, ϵ_surfaces, surfaces, ϵ_charges, charges)
    end
end

Base.show(io::IO, v::PoissonSystem{T, TE}) where{T, TE} = print(io, "PoissonSystem{$T, $TE}, ns: $(ns(v)), nc: $(nc(v))")
Base.copy(p::PoissonSystem{T, TE}) where{T, TE} = PoissonSystem(p.ϵ_medium, copy(p.ϵ_surfaces), copy(p.surfaces), copy(p.ϵ_charges), copy(p.charges))
ns(poisson_sys::PoissonSystem) = length(poisson_sys.surfaces)
nc(poisson_sys::PoissonSystem) = length(poisson_sys.charges)
nt_total(poisson_sys::PoissonSystem) = foldl((x, y) -> x + nt(y), poisson_sys.surfaces, init = 0)

struct PoissonSolver{T, TE, TP}
    A::AbstractArray{TE, 2}
    b::AbstractArray{TE, 1}
    Preconditioner::AbstractArray{TP, 2}
    max_error::T
end