function icosahedron(;T::DataType=Float64)

    m = T(sqrt(50.0 - 10√5) / 10.0)
    n = T(sqrt(50.0 + 10√5) / 10.0)

    @assert m^2 + n^2 ≈ 1.0

    # the vertices of a regular icosahedron
    verts = Vector{Tuple{T, T, T}}()
    for i in [+1, -1]
        for j in [+1, -1]
            push!(verts, (i * m, 0.0, j * n))
            push!(verts, (0.0, j * n, i * m))
            push!(verts, (j * n, i * m, 0.0))
        end
    end

    # the faces of a regular icosahedron
    faces = Vector{Tuple{Int, Int, Int}}()
    r0 = T(1.0514622242382672)
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

    return verts, faces, get_edges(faces)
end

function mesh_sphere(verts::Vector{Tuple{T, T, T}}, faces::Vector{Tuple{Int, Int, Int}}, edges::Vector{Tuple{Int, Int}}) where{T}
    verts2 = copy(verts)

    faces2 = Vector{Tuple{Int, Int, Int}}()

    for face in faces
        a, b, c = face
        v1 = (verts[a] .+ verts[b]) ./ norm(verts[a] .+ verts[b])
        v2 = (verts[b] .+ verts[c]) ./ norm(verts[b] .+ verts[c])
        v3 = (verts[c] .+ verts[a]) ./ norm(verts[c] .+ verts[a])

        d = add_vert!(verts2, v1)
        e = add_vert!(verts2, v2)
        f = add_vert!(verts2, v3)

        push!(faces2, (a, d, f))
        push!(faces2, (b, e, d))
        push!(faces2, (c, f, e))
        push!(faces2, (d, e, f))
    end

    edges2 = get_edges(faces2)

    return verts2, faces2, edges2
end

function sphere(n::Int; T::DataType=Float64)
    verts, faces, edges = icosahedron(T=T)
    for i in 1:n - 1
        verts, faces, edges = mesh_sphere(verts, faces, edges)
    end

    Verteices = Vector{Vertex{T}}()
    for vert in verts
        push!(Verteices, Vertex(vert, vert))
    end

    return Model(Verteices, SVector{3, Int}.(faces))
end