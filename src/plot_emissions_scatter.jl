# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSV of model outputs
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase # get mean

output_dir = "output/default"
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
greenland = DataFrame(CSVFiles.load(joinpath(output_dir, "greenland.csv")))
antarctic = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
gmslr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))


colors = cgrad(:vik100, [0.1, 0.5, 0.9])

idx2100 = findfirst(names(gmslr) .== "2100")
idx2020 = findfirst(names(gmslr) .== "2020")
cum_emissions = [sum(emissions[i, idx2020:idx2100]) for i in 1:nrow(gmslr)]

fig1 = Figure(resolution= (350, 300), fontsize=12, figure_padding=1)
axgmsl = Axis(fig[1,1], xlabel="Global Mean Temperature Anomaly (°C)", ylabel="Global Mean Sea Level Anomaly (m)")


plt =  Makie.scatter!(axgmsl, temperature[:,idx2100  ], gmslr[:, idx2100], color=cum_emissions, colormap=colors, markersize=6)

Colorbar(fig[1,2], plt)

CairoMakie.save("figures/test_plot.png", fig)


idx = findall(cum_emissions .< 1500)

fig2 = Figure(resolution= (350, 300), fontsize=12, figure_padding=15)
axgmsl = Axis(fig2[1,1], xlabel="Global Mean Temperature Anomaly (°C)", ylabel="Global Mean Sea Level Anomaly (m)")
for i in idx
    Makie.lines!(axgmsl, 2020:2100, Vector(gmslr[i, idx2020:idx2100]), color=parameters[i, :t_peak], colorrange=(2030, 2050), colormap=cgrad(:Reds))
end
xlims!(axgmsl, (2020, 2100))
axright = Axis(fig2[1, 2], limits=((0, nothing), (0, 2.0)))
linkyaxes!(axgmsl, axright)
Makie.density!(axright, gmslr[:, idx2100],direction=:y)
hidedecorations!(axright)
hidespines!(axright)

Colorbar(fig2[2, 1:2], limits=(2030, 2050), colormap=:Reds, vertical=:false,    label="Year In Which Emissions Peak", flipaxis=:false)
colsize!(fig2.layout, 2, Auto(0.25))
colgap!(fig2.layout, 1, Relative(0.0))
CairoMakie.save("figures/test_plot.png", fig2)
