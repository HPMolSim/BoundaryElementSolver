using BoundaryElementSolver
using OhMyThreads
using CSV, DataFrames

function main()
    ϵ_0 = 1.0
    ϵ_1 = 10.0

    L = 1.1

    point_charge_1 = PointCharge(0.0, 0.0, L, 1.0)
    point_charge_2 = PointCharge(0.0, 0.0, -L, -1.0)

    bisys = BiSphereSystem(1.0, -1.0, 1.0, 1.0, L, -L, ϵ_0, ϵ_1, ϵ_1)

    CSV.write("data/double_sphere_exact.csv", DataFrame(num_tris = Int[], E_exact = Float64[]))

    for s in 2:14
        surface_1 = sphere_surf(s, center = (0.0, 0.0, L))
        surface_2 = sphere_surf(s, center = (0.0, 0.0, -L))
        poisson_sys0 = PoissonSystem(ϵ_0, [ϵ_1, ϵ_1], [surface_1, surface_2], [ϵ_1, ϵ_1], [point_charge_1, point_charge_2])

        positions = []
        for surf in poisson_sys0.surfaces
            for tri in surf.tris
                push!(positions, tri.r)
            end
        end

        b_phi = zeros(Float64, length(positions))
        @tasks for i in 1:length(positions)
            x, y, z = positions[i]
            b_phi[i] = bisphere_phi(x, y, z, bisys)
        end

        E_exact = energy(poisson_sys0, b_phi)
        num_tris = BoundaryElementSolver.nt_total(poisson_sys0)

        @info "s = $s, E_exact = $E_exact"

        df = DataFrame(num_tris = num_tris, E_exact = E_exact)
        CSV.write("data/double_sphere_exact.csv", df, append = true)
    end
end

main()