using BoundaryElementSolver
using CSV, DataFrames

# consider a single sphere with r = 1.0, with a charge of 1.0 at the center
# ϵ_0 = 80.0 and ϵ_1 = 2.0
# the potential at the surface of the sphere is 1.0 / 4πϵ₀ = 0.0009947183943243459
# the energy of the system is given by (80 - 2) / 2 * 1.0^2 / 4π / 2 / 80.0
function main()

    CSV.write("single_sphere.csv", DataFrame(num_tris = Int[], E = Float64[], errors = Float64[]))

    ϵ_0 = 80.0
    ϵ_1 = 2.0
    E_exact = - (80 - 2) / 2 * 1.0^2 / 4π / 2 / 80.0

    num_tris = Vector{Int}()
    errors = Vector{Float64}()
    # different number of triangles
    E_array = Vector{Float64}()
    for n in 2:7
        point_charge = PointCharge(0.0, 0.0, 0.0, 1.0)
        surf = sphere_surf(n)
        poisson_sys = PoissonSystem(ϵ_0, [ϵ_1], [surf], [ϵ_1], [point_charge])
        φ = solve(poisson_sys)
        E = energy(poisson_sys, φ)
        push!(E_array, E)
        push!(errors, abs(E - E_exact) / abs(E_exact))
        push!(num_tris, length(surf.tris))
        println("num_tris = $(length(surf.tris)), E = $E, error = $(abs(E - E_exact)), relative error = $(abs(E - E_exact) / E_exact)")
        CSV.write("single_sphere.csv", DataFrame(num_tris = num_tris, E = E, errors = errors), append = true)
    end
end

main()