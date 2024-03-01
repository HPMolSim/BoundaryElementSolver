@testset "norm of sphere" begin
    for dt in [Float32, Float64]
        model = load_sphere("n80", data_type = dt)
        obj = Model2Object(model)

        for i in 1:length(obj.tris)
            tri = obj.tris[i]
            n = tri.n
            @test norm(tri.n) ≈ 1.0

            # test the direct of the norms here
        end
    end
end