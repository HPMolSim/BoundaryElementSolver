using BoundaryElementSolver
using GLMakie
using LinearAlgebra

# consider two spheres, one at (0, 0, 2), one at (0, 0, -2), both with r = 1.0
# with center charges of 1.0 and 1.0 respectively

ϵ_0 = 1.0
ϵ_1 = 10.0

bisys = BiSphereSystem(1.0, -1.0, 1.0, 1.0, 1.1, -1.1, ϵ_0, ϵ_1, ϵ_1)
bi_E = bisphere_energy(bisys)

point_charge_1 = PointCharge(0.0, 0.0, 1.1, 10.0)
point_charge_2 = PointCharge(0.0, 0.0, -1.1, -1.0)
s = sphere(4, r = 1.0)
surface_1 = Model2Surface(s, center = (0.0, 0.0, 1.1))
surface_2 = Model2Surface(s, center = (0.0, 0.0, -1.1))
poisson_sys = PoissonSystem(ϵ_0, [ϵ_1, ϵ_1], [surface_1, surface_2], [ϵ_1, ϵ_1], [point_charge_1, point_charge_2])

ϕ = solve(poisson_sys)
surface_phi = surface_potential(poisson_sys, ϕ)

bi_phi = Vector{Float64}()
for (i, surface) in enumerate(poisson_sys.surfaces)
    for (k, tri) in enumerate(surface.tris)
        x, y, z = tri.r
        # c = z > 0 ? 2.0 : -2.0
        # zp = z - c
        # rp = [0, 0, c + 1e-12] .+ [x, y, zp] ./ norm([x, y, zp])
        push!(bi_phi, bisphere_phi(x, y, z, bisys))
    end
end

display_potential(poisson_sys, ϕ; markersize = 10)