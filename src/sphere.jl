function icosahedron(r::T) where{T}

    m = T(sqrt(50.0 - 10√5) / 10.0)
    n = T(sqrt(50.0 + 10√5) / 10.0)

    @assert m^2 + n^2 ≈ 1.0

    # the vertices of a regular icosahedron
    verts = Vector{Tuple{T, T, T}}()
    for i in [+1, -1]
        for j in [+1, -1]
            push!(verts, (r * i * m, zero(T), r * j * n))
            push!(verts, (zero(T), r * j * n, r * i * m))
            push!(verts, (r * j * n, r * i * m, zero(T)))
        end
    end

    # the faces of a regular icosahedron
    faces = Vector{Tuple{Int, Int, Int}}()
    r0 = T(1.0514622242382672) * r
    for i in 1:12
        for j in i + 1:12
            for k in j + 1:12
                r1 = norm(verts[i] .- verts[j])
                r2 = norm(verts[i] .- verts[k])
                r3 = norm(verts[k] .- verts[j])
                if r1 ≈ r2 ≈ r3 ≈ r0 push!(faces, (i, j, k)) end
            end
        end
    end

    return verts, faces
end

function mesh_sphere(verts::Vector{Tuple{T, T, T}}, faces::Vector{Tuple{Int, Int, Int}}, r::T) where{T}
    verts2 = copy(verts)

    faces2 = Vector{Tuple{Int, Int, Int}}()

    for face in faces
        a, b, c = face
        v1 = (verts[a] .+ verts[b]) ./ norm(verts[a] .+ verts[b]) .* r
        v2 = (verts[b] .+ verts[c]) ./ norm(verts[b] .+ verts[c]) .* r
        v3 = (verts[c] .+ verts[a]) ./ norm(verts[c] .+ verts[a]) .* r

        push!(verts2, v1)
        push!(verts2, v2)
        push!(verts2, v3)

        d = length(verts2) - 2
        e = length(verts2) - 1
        f = length(verts2)

        push!(faces2, (a, d, f))
        push!(faces2, (b, e, d))
        push!(faces2, (c, f, e))
        push!(faces2, (d, e, f))
    end

    return verts2, faces2
end

"""
    sphere(n::Int; r::T = 1.0) where{T}

Constructs a sphere model with the specified number of subdivisions.

# Arguments
- `n::Int`: The number of subdivisions to perform on the initial icosahedron, number of faces is `20 * 4^n`.
- `r::T = 1.0`: The radius of the sphere. Default is `1.0`.

# Returns
A `Model` object representing the sphere.

"""
function sphere(n::Int; r::T = 1.0) where{T}
    verts, faces = icosahedron(r)
    for i in 1:n - 1
        verts, faces = mesh_sphere(verts, faces, r)
    end

    Verteices = Vector{Vertex{T}}()
    for vert in verts
        push!(Verteices, Vertex(vert, vert))
    end

    return Model(Verteices, SVector{3, Int}.(faces))
end

function sphere_surf(n::Int; r::T = 1.0, center::NTuple{3, T} = (0.0, 0.0, 0.0)) where{T}
    verts, faces = icosahedron(r)
    for i in 1:n - 1
        verts, faces = mesh_sphere(verts, faces, r)
    end

    n_verts = Vector{Vertex{T}}()
    for vert in verts
        push!(n_verts, Vertex(vert, vert))
    end

    tris = Vector{Triangle{T}}()
    for face in faces
        a, b, c = face
        tri = Triangle(n_verts[a], n_verts[b], n_verts[c])
        tri.r = tri.r ./ norm(tri.r) .* r
        push!(tris, tri)
    end

    surf = MeshedSurface(tris)

    return shift!(surf, center...)
end