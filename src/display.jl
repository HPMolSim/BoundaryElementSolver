function display_model(model::Model{T}) where{T}
    plotly()
    verts = model.verts
    edges = get_edges(model.faces)

    x = [v.p[1] for v in verts]
    y = [v.p[2] for v in verts]
    z = [v.p[3] for v in verts]
    fig = scatter(x, y, z, seriestype = :scatter, label = "")

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