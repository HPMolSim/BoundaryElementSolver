@inline function Coulomb(vert_i::SVector{3, T}, vert_j::SVector{3, T}) where{T}
    return one(T) / norm(vert_i - vert_j) / T(4π)
end

@inline function partial_ni_Coulomb(vert_i::SVector{3, T}, vert_j::SVector{3, T}, n_i::SVector{3, T}) where{T}
    return - dot(vert_i - vert_j, n_i) / norm(vert_i - vert_j)^3 / T(4π)
end

@inline function partial_nj_Coulomb(vert_i::SVector{3, T}, vert_j::SVector{3, T}, n_j::SVector{3, T}) where{T}
    return dot(vert_i - vert_j, n_j) / norm(vert_i - vert_j)^3 / T(4π)
end

@inbounds function energy(poisson_sys::PoissonSystem{T, TE}, ϕ::Vector{TE}) where{T, TE}

    surfaces = poisson_sys.surfaces
    charges = poisson_sys.charges
    ϵ_m = poisson_sys.ϵ_medium
    ϵ_c = poisson_sys.ϵ_charges
    ϵ_s = poisson_sys.ϵ_surfaces

    Nt_total = nt_total(poisson_sys)
    E_pp = zero(TE) # energy between charges and charges
    E_ps = zero(TE) # energy between charges and surfaces

    # for (i, q_i) in enumerate(charges)
    #     for (j, q_j) in enumerate(charges)
    #         (i == j) && continue # skip the self-interaction
    #         E_pp += q_i.q * q_j.q * Coulomb(q_i.r, q_j.r) / ϵ_c[i]
    #     end
    # end

    for (i, q_i) in enumerate(charges)
        jl0 = 0
        for (j, surface_j) in enumerate(surfaces)
            t = tmapreduce(+, 1:nt(surface_j); ntasks = nthreads()) do l
                tri_k = surface_j.tris[l]
                jl = jl0 + l
                γ_ij = TE((ϵ_m - ϵ_s[j]) / ϵ_c[i])
                q_i.q * γ_ij * ϕ[jl] * partial_nj_Coulomb(q_i.r, tri_k.r, tri_k.n) * tri_k.a
            end
            E_ps += t
            jl0 += nt(surface_j)
        end
    end

    @debug "E_pp: $(E_pp / 2), E_ps: $(E_ps / 2)"

    return (E_pp + E_ps) / 2
end

@inbounds function Poisson_A(poisson_sys::PoissonSystem{T, TE}) where{T, TE}

    surfaces = poisson_sys.surfaces
    ϵ_m = poisson_sys.ϵ_medium
    ϵ_s = poisson_sys.ϵ_surfaces

    Nt_total = nt_total(poisson_sys)
    A = zeros(TE, Nt_total, Nt_total) + I

    ik0 = 0

    for (i, surface_i) in enumerate(surfaces)
        jl0 = 0
        for (j, surface_j) in enumerate(surfaces)
            α_ij = TE(2 * (ϵ_m - ϵ_s[j]) / (ϵ_m + ϵ_s[i]))
            for (k, tri_k) in enumerate(surface_i.tris)
                ik = ik0 + k
                @tasks for l in 1:nt(surface_j)
                    if !(i == j && k == l)
                        tri_l = surface_j.tris[l]
                        A[ik, jl0 + l] = - α_ij * partial_nj_Coulomb(tri_k.r, tri_l.r, tri_l.n) * tri_l.a
                    end
                end
            end
            jl0 += nt(surface_j)
        end
        ik0 += nt(surface_i)
    end

    return A
end

@inbounds function Poisson_A_ss(poisson_sys::PoissonSystem{T, TE}) where{T, TE}

    surfaces = poisson_sys.surfaces
    ϵ_m = poisson_sys.ϵ_medium
    ϵ_s = poisson_sys.ϵ_surfaces

    Nt_total = nt_total(poisson_sys)
    A = zeros(TE, Nt_total, Nt_total)
    
    ik0 = 0
    for (i, surface_i) in enumerate(surfaces)
        α_ii = TE(2 * (ϵ_m - ϵ_s[i]) / (ϵ_m + ϵ_s[i]))
        for k in 1:nt(surface_i)
            A[ik0 + k, ik0 + k] = 1 + α_ii / 2
        end
        ik0 += nt(surface_i)
    end

    ik0 = 0
    for (i, surface_i) in enumerate(surfaces)
        jl0 = 0
        for (j, surface_j) in enumerate(surfaces)
            α_ij = TE(2 * (ϵ_m - ϵ_s[j]) / (ϵ_m + ϵ_s[i]))
            @tasks for k in 1:nt(surface_i)
                tri_k = surface_i.tris[k]
                ik = ik0 + k
                for l in 1:nt(surface_j)
                    if !(i == j && k == l)
                        tri_l = surface_j.tris[l]
                        t = - α_ij * partial_nj_Coulomb(tri_k.r, tri_l.r, tri_l.n) * tri_l.a
                        A[ik, jl0 + l] += t
                        (i == j) && (A[ik, ik] -= t)
                    end
                end
            end
            jl0 += nt(surface_j)
        end
        ik0 += nt(surface_i)
    end

    return A
end

@inbounds function Poisson_b(poisson_sys::PoissonSystem{T, TE}) where{T, TE}

    surfaces = poisson_sys.surfaces
    charges = poisson_sys.charges
    ϵ_m = poisson_sys.ϵ_medium
    ϵ_c = poisson_sys.ϵ_charges

    Nt_total = nt_total(poisson_sys) # total number of triangles
    b = zeros(TE, Nt_total)

    ik0 = 0
    for surface in surfaces
        @tasks for k in 1:nt(surface)
            tri = surface.tris[k]
            ik = ik0 + k
            for (j, charge) in enumerate(charges)
                βi  = TE(2 / (ϵ_m + ϵ_c[j]))
                b[ik] += βi * charge.q * Coulomb(tri.r, charge.r)
            end
        end
        ik0 += nt(surface)
    end

    return b
end

function solve(sys::PoissonSystem{T, TE}; ss::Bool=false, kwargs...) where{T, TE}
    A = ss ? Poisson_A_ss(sys) : Poisson_A(sys)
    b = Poisson_b(sys)
    x, _ = gmres(A, b; kwargs...)
    return x
end

function solve_ss(sys::PoissonSystem{T, TE}; kwargs...) where{T, TE}
    A = Poisson_A_ss(sys)
    b = Poisson_b(sys)
    x, _ = gmres(A, b; kwargs...)
    return x
end

function surface_potential(sys::PoissonSystem{T, TE}, ϕ::Vector{T}) where{T, TE}

    surfaces = sys.surfaces
    Nt_total = nt_total(sys)
    potentials = zeros(TE, Nt_total, 4)

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