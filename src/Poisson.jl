@inline function Coulomb(vert_i::SVector{3, T}, vert_j::SVector{3, T}) where{T}
    return one(T) / norm(vert_i - vert_j) / T(4π)
end

@inline function partial_ni_Coulomb(vert_i::SVector{3, T}, vert_j::SVector{3, T}, n_i::SVector{3, T}) where{T}
    return - dot(vert_i - vert_j, n_i) / norm(vert_i - vert_j)^3 / T(4π)
end

@inbounds function energy(poisson_sys::PoissonSystem{T}, ϕ::Vector{T}) where{T}

    surfaces = poisson_sys.surfaces
    charges = poisson_sys.charges
    ϵ_m = poisson_sys.ϵ_medium
    ϵ_c = poisson_sys.ϵ_charges
    ϵ_s = poisson_sys.ϵ_surfaces

    Nt_total = nt_total(poisson_sys)
    E = zero(T) # total energy

    for (i, q_i) in enumerate(charges)
        for (j, q_j) in enumerate(charges)
            (i == j) && continue # skip the self-interaction
            E += q_i.q * q_j.q * Coulomb(q_i.r, q_j.r) / ϵ_c[i]
        end
    end

    for (i, q_i) in enumerate(charges)
        jl = 0
        for (j, surface_j) in enumerate(surfaces)
            for (l, tri_k) in enumerate(surface_j.tris)
                jl += 1
                γ_ij = T((ϵ_m - ϵ_s[j]) / ϵ_c[i])
                E -= q_i.q * γ_ij * ϕ[jl] * partial_ni_Coulomb(tri_k.r, q_i.r, tri_k.n) * tri_k.a
            end
        end
    end


    return E / 2
end

@inbounds function Poisson_A(poisson_sys::PoissonSystem{T}) where{T}

    surfaces = poisson_sys.surfaces
    ϵ_m = poisson_sys.ϵ_medium
    ϵ_s = poisson_sys.ϵ_surfaces

    Nt_total = nt_total(poisson_sys)
    A = zeros(T, Nt_total, Nt_total) + I

    ik0 = 0

    for (i, surface_i) in enumerate(surfaces)
        jl0 = 0
        for (j, surface_j) in enumerate(surfaces)
            α_ij = T(2 * (ϵ_m - ϵ_s[j]) / (ϵ_m + ϵ_s[i]))
            for (k, tri_k) in enumerate(surface_i.tris)
                ik = ik0 + k
                for (l, tri_l) in enumerate(surface_j.tris)
                    (i == j && k == l) && continue # skip the self-interaction
                    jl = jl0 + l
                    A[ik, jl] += - α_ij * partial_ni_Coulomb(tri_l.r, tri_k.r, tri_l.n) * tri_l.a
                end
            end
            jl0 += nt(surface_j)
        end
        ik0 += nt(surface_i)
    end

    return A
end

@inbounds function Poisson_b(poisson_sys::PoissonSystem{T}) where{T}

    surfaces = poisson_sys.surfaces
    charges = poisson_sys.charges
    ϵ_m = poisson_sys.ϵ_medium
    ϵ_c = poisson_sys.ϵ_charges

    Nt_total = nt_total(poisson_sys) # total number of triangles
    b = zeros(T, Nt_total)

    ik = 0
    for (i, surface) in enumerate(surfaces)
        for (k, tri) in enumerate(surface.tris)
            ik += 1
            for (j, charge) in enumerate(charges)
                βi  = T(2 / (ϵ_m + ϵ_c[j]))
                b[ik] += βi * charge.q * Coulomb(tri.r, charge.r)
            end
        end
    end

    return b
end

function solve(sys::PoissonSystem{T}; kwargs...) where{T}
    A = Poisson_A(sys)
    b = Poisson_b(sys)
    x, _ = gmres(A, b; kwargs...)
    return x
end

function surface_potential(sys::PoissonSystem{T}, ϕ::Vector{T}) where{T}

    surfaces = sys.surfaces
    Nt_total = nt_total(sys)
    potentials = zeros(T, Nt_total, 4)

    ik = 0
    for (i, surface) in enumerate(surfaces)
        for (k, tri) in enumerate(surface.tris)
            ik += 1
            x, y, z = tri.r
            potentials[ik, 1] = x
            potentials[ik, 2] = y
            potentials[ik, 3] = z
            potentials[ik, 4] = ϕ[ik]
        end
    end

    return potentials
end