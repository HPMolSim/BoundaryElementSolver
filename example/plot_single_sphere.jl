using CairoMakie, LaTeXStrings
using CSV, DataFrames

df = CSV.read("single_sphere.csv", DataFrame)

E_exact = - (80 - 2) / 2 * 1.0^2 / 4π / 2 / 80.0
num_tris = df.num_tris

E_0 = df.E
E_1 = [2 * E_0[i + 1] - E_0[i] for i in 1:length(E_0) - 1]
E_2 = [(4 * E_1[i + 1] - E_1[i]) / 3 for i in 1:length(E_1) - 1]

ms = 15
fs = 20
marker_style = [:circle, :diamond, :utriangle, :dtriangle]
colors = [:blue, :green, :red, :purple, :brown]

begin
    fig = Figure()
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10, title = "single sphere", xlabel = "number of triangles", ylabel = "relative error")

    scatterlines!(ax, num_tris, abs.(E_0 .- E_exact) ./ abs(E_exact), label = "CC", markersize = ms, linestyle = :dash, marker = marker_style[1], color = colors[1])
    scatterlines!(ax, num_tris[2:end], abs.(E_1 .- E_exact) ./ abs(E_exact), label = "CC-1", markersize = ms, linestyle = :dash, marker = marker_style[2], color = colors[2])
    scatterlines!(ax, num_tris[3:end], abs.(E_2 .- E_exact) ./ abs(E_exact), label = "CC-2", markersize = ms, linestyle = :dash, marker = marker_style[3], color = colors[3])


    lines!(ax, num_tris[end - 1:end], num_tris[end - 1:end].^(-0.5) ./ 1, color = :black)
    lines!(ax, num_tris[end - 1:end], num_tris[end - 1:end].^(-1) * 10, color = :black)
    lines!(ax, num_tris[end - 1:end], num_tris[end - 1:end].^(-1.5) * 50, color = :black)

    text!(ax, (1e5, 1e-2), text = L"N^{-0.5}", color = :black, align = (:right, :top), fontsize = fs)
    text!(ax, (1e5, 10^(-3.3)), text = L"N^{-1}", color = :black, align = (:right, :top), fontsize = fs)
    text!(ax, (1e5, 1e-5), text = L"N^{-1.5}", color = :black, align = (:right, :top), fontsize = fs)

    axislegend(ax, position = :lb)
end
fig

save("single_sphere.png", fig)
