function display_model(model::Model{T}) where{T}
    plotly()
    verts = model.verts
    edges = get_edges(model.faces)

    x = [v.p[1] for v in verts]
    y = [v.p[2] for v in verts]
    z = [v.p[3] for v in verts]
    fig = scatter(x, y, z, seriestype = :scatter, label = "", color = :blue)

    for edge in edges
        i, j = edge
        plot!(fig, 
            [verts[i].p[1], verts[j].p[1]],
            [verts[i].p[2], verts[j].p[2]],
            [verts[i].p[3], verts[j].p[3]],
            label = "", color = :black)
    end

    return fig
end

function display_surface!(fig, surface::MeshedSurface{T}) where{T}
    tris = surface.tris
    x = [tri.r[1] for tri in tris]
    y = [tri.r[2] for tri in tris]
    z = [tri.r[3] for tri in tris]

    return scatter!(fig, x, y, z, seriestype = :scatter, label = "", color = :blue, markersize = 2)
end

function display_surface(surface::MeshedSurface{T}) where{T}
    plotly()
    fig = plot()
    display_surface!(fig, surface)
    return fig
end

function display_system(poisson_sys::PoissonSystem{T}) where{T}
    plotly()
    fig = plot()
    for surface in poisson_sys.surfaces
        display_surface!(fig, surface)
    end

    for charge in poisson_sys.charges
        scatter!(fig, [charge.r[1]], [charge.r[2]], [charge.r[3]], seriestype = :scatter, label = "", color = :red, markersize = 2)
    end

    return fig
end