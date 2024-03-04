function display_model(model::Model{T}) where{T}
    plotly()
    verts = model.verts
    edges = get_edges(model)

    x = [v[1][1] for v in verts]
    y = [v[1][2] for v in verts]
    z = [v[1][3] for v in verts]
    fig = scatter3d(x, y, z, seriestype = :scatter, label = "")

    for edge in edges
        i, j = edge
        plot!(fig, 
            [verts[i][1][1], verts[j][1][1]],
            [verts[i][1][2], verts[j][1][2]],
            [verts[i][1][3], verts[j][1][3]],
            label = "", color = :black)
    end

    return fig
end