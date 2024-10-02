using BoundaryElementSolver
using CSV, DataFrames
using LinearAlgebra

function main()
    # consider two spheres, one at (0, 0, L), one at (0, 0, -L), both with r = 1.0
    # with center charges of 1.0 and 1.0 respectively

    ϵ_0 = 1.0
    ϵ_1 = 10.0

    L = 1.1

    point_charge_1 = PointCharge(0.0, 0.0, L, 1.0)
    point_charge_2 = PointCharge(0.0, 0.0, -L, -1.0)

    bisys = BiSphereSystem(1.0, -1.0, 1.0, 1.0, L, -L, ϵ_0, ϵ_1, ϵ_1)
    surface_1 = sphere_surf(8, center = (0.0, 0.0, L))
    surface_2 = sphere_surf(8, center = (0.0, 0.0, -L))
    poisson_sys0 = PoissonSystem(ϵ_0, [ϵ_1, ϵ_1], [surface_1, surface_2], [ϵ_1, ϵ_1], [point_charge_1, point_charge_2])

    positions = []
    for surf in poisson_sys0.surfaces
        for tri in surf.tris
            push!(positions, tri.r)
        end
    end

    b_phi = Float64[]
    for pos in positions
        x, y, z = pos
        push!(b_phi, bisphere_phi(x, y, z, bisys))
    end

    E_exact = energy(poisson_sys0, b_phi)

    CSV.write("data/double_sphere.csv", DataFrame(num_tris = Int[], E_exact = Float64[], E_cc = Float64[], E_ss = Float64[]))

    for s in 2:6
        surface_1 = sphere_surf(s, center = (0.0, 0.0, L))
        surface_2 = sphere_surf(s, center = (0.0, 0.0, -L))
        poisson_sys = PoissonSystem(ϵ_0, [ϵ_1, ϵ_1], [surface_1, surface_2], [ϵ_1, ϵ_1], [point_charge_1, point_charge_2])

        ϕ = solve(poisson_sys)
        E = energy(poisson_sys, ϕ)

        ϕ_ss = solve_ss(poisson_sys)
        E_ss = energy(poisson_sys, ϕ_ss)

        @show E, E_exact, abs(E - E_exact) / abs(E_exact), E_ss, abs(E_ss - E_exact) / abs(E_exact)

        ntt = BoundaryElementSolver.nt_total(poisson_sys)

        CSV.write("data/double_sphere.csv", DataFrame(num_tris = ntt, E_exact = E_exact, E_cc = E, E_ss = E_ss), append = true)
    end
end

main()