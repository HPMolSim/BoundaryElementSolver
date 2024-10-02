using BoundaryElementSolver
using CSV, DataFrames
using LinearAlgebra

using CairoMakie, LaTeXStrings

function main()
    # consider two spheres, one at (0, 0, L), one at (0, 0, -L), both with r = 1.0
    # with center charges of 1.0 and 1.0 respectively

    ϵ_0 = 1.0
    ϵ_1 = 10.0

    L = 1.1

    point_charge_1 = PointCharge(0.0, 0.0, L, 1.0)
    point_charge_2 = PointCharge(0.0, 0.0, -L, -1.0)

    bisys = BiSphereSystem(1.0, -1.0, 1.0, 1.0, L, -L, ϵ_0, ϵ_1, ϵ_1)
    surface_1 = sphere_surf(9, center = (0.0, 0.0, L))
    surface_2 = sphere_surf(9, center = (0.0, 0.0, -L))
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

# E_cc_ss_1 = [(4 * E_ss_array[i + 1] - E_ss_array[i]) / 3 for i in 1:length(E_ss_array) - 1]
# E_cc_ss_2 = [(16 * E_cc_ss_1[i + 1] - E_cc_ss_1[i]) / 15 for i in 1:length(E_cc_ss_1) - 1]

# begin
#     num_tris = [2 * 20 * 4^i for i in 2:5]
#     errors_cc = abs.(E_array .- E_exact) ./ abs.(E_exact)
#     errors_ss = abs.(E_ss_array .- E_exact) ./ abs.(E_exact)
#     errors_cc_ss_1 = abs.(E_cc_ss_1 .- E_exact) ./ abs.(E_exact)
#     errors_cc_ss_2 = abs.(E_cc_ss_2 .- E_exact) ./ abs.(E_exact)

#     fig = Figure()
#     ax = Axis(fig[1, 1], xscale = log10, yscale = log10, title = "double sphere", xlabel = "number of triangles", ylabel = "relative error")
#     scatterlines!(ax, num_tris, errors_cc, label = "CC", markersize = 15, linestyle = :dash)
#     scatterlines!(ax, num_tris, errors_ss, label = "SS", markersize = 15, linestyle = :dash)
#     scatterlines!(ax, num_tris[2:end], errors_cc_ss_1, label = "CC-SS-1", markersize = 15, linestyle = :dash)
#     scatterlines!(ax, num_tris[3:end], errors_cc_ss_2, label = "CC-SS-2", markersize = 15, linestyle = :dash)
#     # lines!(ax, num_tris, num_tris.^(-0.5) .* 13, color = :red, label = L"N^{-0.5}")
#     axislegend(ax, position = :lb)
# end
# fig

# save("figs/double_sphere.png", fig)