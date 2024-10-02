using CSV, DataFrames
using CairoMakie, LaTeXStrings

df = CSV.read("data/double_sphere.csv", DataFrame)
num_tris = df.num_tris
# E_exact = df.E_exact[1]
E_exact = 0.03470251241072653
E_cc = df.E_cc
E_ss = df.E_ss

E_cc_ss_1 = [(4 * E_ss[i + 1] - E_ss[i]) / 3 for i in 1:length(E_ss) - 1]
E_cc_ss_2 = [(16 * E_cc_ss_1[i + 1] - E_cc_ss_1[i]) / 15 for i in 1:length(E_cc_ss_1) - 1]

errors_cc = abs.(E_cc .- E_exact) ./ abs.(E_exact)
errors_ss = abs.(E_ss .- E_exact) ./ abs.(E_exact)
errors_cc_ss_1 = abs.(E_cc_ss_1 .- E_exact) ./ abs.(E_exact)
errors_cc_ss_2 = abs.(E_cc_ss_2 .- E_exact) ./ abs.(E_exact)

begin

    fig = Figure()
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10, title = "double sphere", xlabel = "number of triangles", ylabel = "relative error")
    scatterlines!(ax, num_tris, errors_cc, label = "CC", markersize = 15, linestyle = :dash)
    scatterlines!(ax, num_tris, errors_ss, label = "SS", markersize = 15, linestyle = :dash)
    scatterlines!(ax, num_tris[2:end], errors_cc_ss_1, label = "CC-SS-1", markersize = 15, linestyle = :dash)
    scatterlines!(ax, num_tris[3:end], errors_cc_ss_2, label = "CC-SS-2", markersize = 15, linestyle = :dash)
    
    ns = num_tris[end - 1:end]
    lines!(ax, ns, ns.^(-0.5) * errors_cc[end] / ns[2]^(-0.5) * 2 , color = :black)
    lines!(ax, ns, ns.^(-1) * errors_ss[end] / ns[2]^(-1) * 2 , color = :black)
    lines!(ax, ns, ns.^(-2) * errors_cc_ss_1[end] / ns[2]^(-2) * 2 , color = :black)
    lines!(ax, ns, ns.^(-3) * errors_cc_ss_2[end] / ns[2]^(-3) / 2 , color = :black)

    axislegend(ax, position = :lb)
end
fig

save("figs/double_sphere.png", fig)