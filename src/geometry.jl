# the operations about triangular meshes

@inline function tri_area(a::AbstractArray{T, 1}, b::AbstractArray{T, 1}, c::AbstractArray{T, 1}) where T
    normal = cross(b - a, c - a)
    return T(norm(normal) / 2)
end

function Model2Surface(model::Model{T}; center::NTuple{3, T} = (zero(T), zero(T), zero(T))) where{T}
    verts = model.verts
    faces = model.faces
    tris = Vector{Triangle{T}}()
    for face in faces
        a, b, c = face
        push!(tris, Triangle(verts[a], verts[b], verts[c]))
    end

    surf = MeshedSurface(tris)

    return shift!(surf, center...)
end

function get_edges(faces::Vector)
    edges = Vector{Tuple{Int, Int}}()
    for face in faces
        a, b, c = face
        for edge in [(a, b), (b, c), (c, a)]
            if !(edge in edges) && !((edge[2], edge[1]) in edges)
                push!(edges, edge)
            end
        end
    end
    return edges
end

function add_vert!(verts, new_vert)
    if new_vert in verts
        return findfirst(isequal(new_vert), verts)
    else
        push!(verts, new_vert)
        return length(verts)
    end
end

function shift!(m::MeshedSurface{T}, x::T, y::T, z::T) where{T}
    for tri in m.tris
        tri.r = tri.r .+ SVector{3, T}(x, y, z)
    end
    return m
end