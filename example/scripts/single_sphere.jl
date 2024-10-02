using BoundaryElementSolver
using CSV, DataFrames

# consider a single sphere with r = 1.0, with a charge of 1.0 at the center
# ϵ_0 = 80.0 and ϵ_1 = 2.0
# the potential at the surface of the sphere is 1.0 / 4πϵ₀ = 0.0009947183943243459
# the energy of the system is given by (80 - 2) / 2 * 1.0^2 / 4π / 2 / 80.0
function main()

    CSV.write("data/single_sphere.csv", DataFrame(num_tris = Int[], E_exact = Float64[], E_cc = Float64[], E_ss = Float64[]))

    ϵ_0 = 80.0
    ϵ_1 = 2.0
    E_exact = - (80 - 2) / 2 * 1.0^2 / 4π / 2 / 80.0

    num_tris = Vector{Int}()
    errors = Vector{Float64}()
    # different number of triangles
    for n in 2:7
        point_charge = PointCharge(0.0, 0.0, 0.0, 1.0)
        surf = sphere_surf(n)
        poisson_sys = PoissonSystem(ϵ_0, [ϵ_1], [surf], [ϵ_1], [point_charge])
        ϕ_cc = solve(poisson_sys)
        E_cc = energy(poisson_sys, ϕ_cc)

        ϕ_ss = solve_ss(poisson_sys)
        E_ss = energy(poisson_sys, ϕ_ss)

        println("num_tris = $(length(surf.tris)), E_cc = $E_cc, relative error cc = $(abs(E_cc - E_exact) / abs(E_exact)), E_ss = $(E_ss), relative error ss = $(abs(E_ss - E_exact) / abs(E_exact))")

        df = DataFrame(num_tris = length(surf.tris), E_exact = E_exact, E_cc = E_cc, E_ss = E_ss)
        CSV.write("data/single_sphere.csv", df, append = true)
    end
end

main()