using BoundaryElementSolver
using Test, LinearAlgebra
using BoundaryElementSolver: Vertex, mid, ns, nc, nt_total

@testset "Vertex" begin
    v1 = Vertex((1.0, 2.0, 3.0), (1.0, 1.0, 1.0))
    v2 = Vertex((1.0, 2.0, 3.0), (2.0, 2.0, 2.0))
    v3 = copy(v1)

    @test v1 == v3
    @test v1 == v2

    @test mid(v1, v2) == v1
    @test mid(v1, v2, v3) == v1
end

@testset "Triangle" begin
    v1 = Vertex((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    v2 = Vertex((1.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    v3 = Vertex((0.0, 1.0, 0.0), (1.0, 1.0, 1.0))

    tri = Triangle(v1, v2, v3)
    @test tri.r ≈ SVector{3, Float64}((1/3, 1/3, 0.0))
    @test tri.n ≈ SVector{3, Float64}((1/sqrt(3), 1/sqrt(3), 1/sqrt(3)))
    @test tri.a ≈ 0.5
end

@testset "PointCharge" begin
    @test PointCharge(1.0, 2.0, 3.0, 1.0) == PointCharge(SVector{3, Float64}((1.0, 2.0, 3.0)), 1.0)
end

@testset "PoissonSystem" begin
    surf_1 = Model2Surface(sphere(2, r=1.0))
    surf_2 = Model2Surface(sphere(3, r=2.0), center = (5.0, 2.0, 3.0))
    pc = PointCharge(5.0, 2.0, 3.0, 1.0)

    sys = PoissonSystem(1.0, [0.0, 0.0], [surf_1, surf_2], [0.0], [pc])
    @test ns(sys) == 2
    @test nc(sys) == 1
    @test nt_total(sys) == 400
end