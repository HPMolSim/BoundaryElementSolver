# the operations about triangular meshes

@inline function tri_area(a::AbstractArray{T, 1}, b::AbstractArray{T, 1}, c::AbstractArray{T, 1}) where T
    normal = cross(b - a, c - a)
    return T(norm(normal) / 2)
end

function Model2Object(model::Model{T}; center::NTuple{3, T} = (zero(T), zero(T), zero(T)), angle::NTuple{3, T} = (zero(T), zero(T), zero(T)), scale::T = one(T), ϵ::T = one(T)) where{T}
    verts = model.verts
    faces = model.faces
    tris = Vector{Triangle{T}}()
    for face in faces
        a, b, c = face
        r = (verts[a].p + verts[b].p + verts[c].p) .* scale / 3
        n = (verts[a].n + verts[b].n + verts[c].n) ./ 3
        n = n ./ (norm(n))
        a = tri_area(verts[a].p, verts[b].p, verts[c].p)
        push!(tris, Triangle(r, n, a))
    end

    center = SVector{3, T}(center)
    angle = SVector{3, T}(angle)

    return MeshedObject(tris, center, angle, ϵ)
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