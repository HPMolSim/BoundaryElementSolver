@testset "single sphere" begin
    ϵ_0 = 80.0
    ϵ_1 = 2.0
    E_exact = - (80 - 2) / 2 * 1.0^2 / 4π / 2 / 80.0

    point_charge = PointCharge(0.0, 0.0, 0.0, 1.0)
    s = sphere(4, r = 1.0)
    surf = Model2Surface(s)
    poisson_sys = PoissonSystem(ϵ_0, [ϵ_1], [surf], [ϵ_1], [point_charge])

    ϕ = solve(poisson_sys)
    E = energy(poisson_sys, ϕ)

    ϕ_ss = solve_ss(poisson_sys)
    E_ss = energy(poisson_sys, ϕ_ss)

    @test isapprox(E, E_exact, atol = 1e-3)
    @test isapprox(E_ss, E_exact, atol = 1e-3)
end