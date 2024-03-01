function convert_vert(vert_lines::Vector{String}; data_type::DataType = Float64)
    verts = Vector{Tuple{Int64, Tuple{SVector{3, data_type}, SVector{3, data_type}}}}()
    for vert_line in vert_lines
        id, nx, ny, nz, x, y, z, u1, u2 = split(vert_line)
        id = parse(Int64, id)
        nx, ny, nz, x, y, z = map(x -> parse(data_type, x), [nx, ny, nz, x, y, z])
        nx, ny, nz = [nx, ny, nz] ./ sqrt(nx^2 + ny^2 + nz^2)
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

function model_from_artifact(artifact_name::String, name::String; data_type::DataType = Float64)
    path_artifact = @artifact_str artifact_name
    path_vert = joinpath(path_artifact, name * "vert.dat")
    path_face = joinpath(path_artifact, name * "face.dat")

    verts, faces = load_model(path_vert, path_face; data_type)

    return Model{data_type}(verts, faces)
end

function load_sphere(name::String; data_type::DataType = Float64)
    return model_from_artifact("geodesic_grid", name; data_type)
end