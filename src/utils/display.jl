function display_model(model::Model{T}) where{T}

    fig = Figure()
    scene = LScene(fig[1, 1])

    verts = model.verts
    edges = get_edges(model.faces)

    x = [v.p[1] for v in verts]
    y = [v.p[2] for v in verts]
    z = [v.p[3] for v in verts]
    scatter!(scene, x, y, z, label = "", color = :blue, markersize = 15)

    for edge in edges
        i, j = edge
        lines!(scene, 
            [verts[i].p[1], verts[j].p[1]],
            [verts[i].p[2], verts[j].p[2]],
            [verts[i].p[3], verts[j].p[3]],
            label = "", color = :black)
    end

    return fig
end

function display_surface!(scene, surface::MeshedSurface{T}, surface_color, markersize) where{T}
    tris = surface.tris
    x = [tri.r[1] for tri in tris]
    y = [tri.r[2] for tri in tris]
    z = [tri.r[3] for tri in tris]

    return scatter!(scene, x, y, z, color = surface_color, markersize = markersize)
end

function display_surface(surface::MeshedSurface{T}; surface_color = :blue, markersize = 10) where{T}
    fig = Figure()
    scene = LScene(fig[1, 1])
    display_surface!(scene, surface, surface_color, markersize)
    return fig
end

function display_system(poisson_sys::PoissonSystem{T}; 
    surface_color = :blue, 
    charge_color = :red, 
    markersize = 10) where{T}
    fig = Figure()
    scene = LScene(fig[1, 1])
    for surface in poisson_sys.surfaces
        display_surface!(scene, surface, surface_color, markersize)
    end

    for charge in poisson_sys.charges
        scatter!(scene, [charge.r[1]], [charge.r[2]], [charge.r[3]], label = "", color = charge_color, markersize = markersize)
    end

    return fig
end

function display_potential(poisson_sys::PoissonSystem{T}, ϕ::Vector{T}; markersize = 10, colormap = :bluesreds) where{T}
    fig = Figure()
    scene = LScene(fig[1, 1])
    ϕp = ϕ ./ maximum(abs.(ϕ))
    potentials = surface_potential(poisson_sys, ϕp)
    x = potentials[:, 1]
    y = potentials[:, 2]
    z = potentials[:, 3]
    c = potentials[:, 4]
    scatter!(scene, x, y, z, color = c, markersize = markersize, colormap = colormap, colorrange = (-1, 1))
    return fig
end