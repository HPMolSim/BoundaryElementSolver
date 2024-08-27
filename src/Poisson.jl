function Coulomb(vert_i::SVector{3, T}, vert_j::SVector{3, T}) where{T}
    return one(T) / norm(vert_i - vert_j)
end

function partial_ni_Coulomb(vert_i::SVector{3, T}, vert_j::SVector{3, T}, n_i::SVector{3, T}) where{T}
    return - dot(vert_i - vert_j, n_i) / norm(vert_i - vert_j)^3
end

function Poisson_A(poisson_sys::PoissonSystem{T}) where{T}

    surfaces = poisson_sys.surfaces
    charges = poisson_sys.charges
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
                for (l, tri_l) in enumerate(surface_j.tris)
                    (i == j && k == l) && continue # skip the self-interaction
                    ik = ik0 + k
                    jl = jl0 + l
                    A[ik, jl] += α_ij * partial_ni_Coulomb(tri_k.r, tri_l.r, tri_k.n) * tri_k.a
                end
            end
            jl0 += nt(surface_j)
        end
        ik0 += nt(surface_i)
    end

    return A
end

function Poisson_b(poisson_sys::PoissonSystem{T}) where{T}

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