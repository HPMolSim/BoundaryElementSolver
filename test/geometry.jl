@testset "norm of sphere" begin
    for T in [Float32, Float64]
        model = sphere(2, r = T(1.0))
        surf = Model2Surface(model)

        for tri in surf.tris
            @test norm(tri.n) ≈ 1.0
            # test the direct of the norms here
        end
    end
end