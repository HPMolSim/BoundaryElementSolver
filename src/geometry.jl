# the operations about triangular meshes

@inline function tri_area(a::AbstractArray{T, 1}, b::AbstractArray{T, 1}, c::AbstractArray{T, 1}) where T
    normal = cross(b - a, c - a)
    return norm(normal) / 2
end

function Model2Object(model::Model{T}; center::NTuple{3, T} = (zero(T), zero(T), zero(T)), angle::NTuple{3, T} = (zero(T), zero(T), zero(T)), scale::T = one(T), ϵ::T = one(T)) where{T}
    verts = model.verts
    faces = model.faces
    tris = Vector{Triangle{T}}()
    for face in faces
        a, b, c = face
        r = (verts[a][1] + verts[b][1] + verts[c][1]) .* scale / 3
        n = (verts[a][2] + verts[b][2] + verts[c][2]) ./ 3
        n = n ./ (norm(n))
        a = tri_area(verts[a][1], verts[b][1], verts[c][1])
        push!(tris, Triangle(r, n, a))
    end

    center = SVector{3, T}(center)
    angle = SVector{3, T}(angle)

    return MeshedObject(tris, center, angle, ϵ)
end