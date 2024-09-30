using BoundaryElementSolver
using GLMakie, LaTeXStrings
using LinearAlgebra

# consider two spheres, one at (0, 0, 2), one at (0, 0, -2), both with r = 1.0
# with center charges of 1.0 and 1.0 respectively

ϵ_0 = 1.0
ϵ_1 = 10.0

L = 1.5

bisys = BiSphereSystem(1.0, -1.0, 1.0, 1.0, L, -L, ϵ_0, ϵ_1, ϵ_1)
bi_E = bisphere_energy(bisys)

point_charge_1 = PointCharge(0.0, 0.0, L, 1.0)
point_charge_2 = PointCharge(0.0, 0.0, -L, -1.0)

E_array = []

for s in 1:5
    surface_1 = sphere_surf(s, center = (0.0, 0.0, L))
    surface_2 = sphere_surf(s, center = (0.0, 0.0, -L))
    poisson_sys = PoissonSystem(ϵ_0, [ϵ_1, ϵ_1], [surface_1, surface_2], [ϵ_1, ϵ_1], [point_charge_1, point_charge_2])

    # A = BoundaryElementSolver.Poisson_A(poisson_sys)
    # b = BoundaryElementSolver.Poisson_b(poisson_sys)

    ϕ = solve(poisson_sys)
    surface_phi = surface_potential(poisson_sys, ϕ)

    E = energy(poisson_sys, ϕ)

    push!(E_array, E)
end

surface_1 = sphere_surf(6, center = (0.0, 0.0, L))
surface_2 = sphere_surf(6, center = (0.0, 0.0, -L))
poisson_sys = PoissonSystem(ϵ_0, [ϵ_1, ϵ_1], [surface_1, surface_2], [ϵ_1, ϵ_1], [point_charge_1, point_charge_2])

bi_phi = Vector{Float64}()
for surf in poisson_sys.surfaces
    for tri in surf.tris
        x, y, z = tri.r
        push!(bi_phi, bisphere_phi(x, y, z, bisys))
    end
end

E_exact = energy(poisson_sys, bi_phi)

begin
    num_tris = [2 * 20 * 4^i for i in 1:5]
    errors = abs.(E_array .- E_exact) ./ E_exact

    fig = Figure()
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10, title = "double sphere", xlabel = "number of triangles", ylabel = "relative error")
    scatterlines!(ax, num_tris, errors)
    lines!(ax, num_tris, num_tris.^(-0.5) .* 13, color = :red, label = L"N^{-0.5}")
    axislegend(ax)
end
fig

save("example/double_sphere.png", fig)