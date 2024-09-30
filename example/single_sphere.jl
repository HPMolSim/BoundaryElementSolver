using BoundaryElementSolver
using GLMakie, LaTeXStrings

# consider a single sphere with r = 1.0, with a charge of 1.0 at the center
# ϵ_0 = 80.0 and ϵ_1 = 2.0
# the potential at the surface of the sphere is 1.0 / 4πϵ₀ = 0.0009947183943243459
# the energy of the system is given by (80 - 2) / 2 * 1.0^2 / 4π / 2 / 80.0

ϵ_0 = 80.0
ϵ_1 = 2.0
E_exact = - (80 - 2) / 2 * 1.0^2 / 4π / 2 / 80.0

point_charge = PointCharge(0.0, 0.0, 0.0, 1.0)
s = sphere(4, r = 1.0)
surf = Model2Surface(s)
poisson_sys = PoissonSystem(ϵ_0, [ϵ_1], [surf], [ϵ_1], [point_charge])

ϕ = solve(poisson_sys)
E = energy(poisson_sys, ϕ)

num_tris = Vector{Int}()
errors = Vector{Float64}()
# different number of triangles
for n in 1:6
    surf = sphere_surf(n)
    poisson_sys = PoissonSystem(ϵ_0, [ϵ_1], [surf], [ϵ_1], [point_charge])
    φ = solve(poisson_sys)
    E = energy(poisson_sys, φ)
    push!(errors, abs(E - E_exact) / abs(E_exact))
    push!(num_tris, length(surf.tris))
    println("num_tris = $(length(surf.tris)), E = $E, error = $(abs(E - E_exact)), relative error = $(abs(E - E_exact) / E_exact)")
end

begin
    fig = Figure()
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10, title = "single sphere", xlabel = "number of triangles", ylabel = "relative error")
    scatterlines!(ax, num_tris, errors)
    lines!(ax, num_tris, num_tris.^(-0.5) ./ 2, color = :red, label = L"N^{-0.5}")
    axislegend(ax)
end
fig

save("example/single_sphere.png", fig)