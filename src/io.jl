function convert_vert(vert_lines::Vector{String}; data_type::DataType = Float64)
    verts = Vector{Tuple{Int64, Tuple{SVector{3, data_type}, SVector{3, data_type}}}}()
    for vert_line in vert_lines
        id, x, y, z, nx, ny, nz, u1, u2 = split(vert_line)
        id = parse(Int64, id)
        x, y, z, nx, ny, nz = map(x -> parse(data_type, x), [x, y, z, nx, ny, nz])
        push!(verts, (id, (SVector{3, data_type}(x, y, z), SVector{3, data_type}(nx, ny, nz))))
    end
    return Dict(verts)
end

function convert_face(face_lines::Vector{String})
    faces = Vector{SVector{3, Int}}()
    for face_line in face_lines
        f1, f2, f3 = split(face_line)
        f1, f2, f3 = map(x -> parse(Int64, x), [f1, f2, f3])
        push!(faces, SVector{3, Int64}(f1, f2, f3))
    end
    return faces
end

function load_model(path_vert::String, path_face::String; data_type::DataType = Float64)

    vert_lines = try
        open(path_vert) do file
            readlines(file)
        end
    catch
        @error "vert file not found."
    end

    face_lines = try
        open(path_face) do file
            readlines(file)
        end
    catch
        @error "face file not found."
    end

    vert_dat = convert_vert(vert_lines; data_type)
    face_dat = convert_face(face_lines)

    return vert_dat, face_dat
end

function model_from_artifact(name::String; data_type::DataType = Float64)
    path_vert = joinpath(artifact"geodesic_grid", name * "vert.dat")
    path_face = joinpath(artifact"geodesic_grid", name * "face.dat")

    verts, faces = load_model(path_vert, path_face; data_type)

    return Model{data_type}(verts, faces)
end