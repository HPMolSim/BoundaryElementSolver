using BoundaryElementSolver
using BoundaryElementSolver: icosahedron, mesh_sphere
using Test

@testset "sphere" begin
    for T in [Float32, Float64]
        for r in [1.0, 2.0]
            verts, faces = icosahedron(T(r))
            @test length(verts) == 12
            @test length(faces) == 20
            for v in verts
                @test sqrt(v[1]^2 + v[2]^2 + v[3]^2) ≈ r
            end
        end
    end

    for T in [Float32, Float64]
        for r in [1.0, 2.0]
            verts, faces = icosahedron(T(r))
            verts2, faces2 = mesh_sphere(verts, faces, T(r))
            @test length(verts2) == 42
            @test length(faces2) == 80
            for v in verts2
                @test sqrt(v[1]^2 + v[2]^2 + v[3]^2) ≈ r
            end

            model = sphere(2, r=T(r))
            for v in model.verts
                @test sqrt(v.p[1]^2 + v.p[2]^2 + v.p[3]^2) ≈ r
                @test v.p ./ T(r) ≈ v.n
            end
        end
    end
end