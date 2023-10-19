# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSVs
using XLSX # read XL files
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase

output_dir = "output/default"
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
gmslr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
antarctic = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))
gsic = DataFrame(CSVFiles.load(joinpath(output_dir, "gsic.csv")))
greenland = DataFrame(CSVFiles.load(joinpath(output_dir, "greenland.csv")))
lw_storage = DataFrame(CSVFiles.load(joinpath(output_dir, "lw_storage.csv")))
thermal_expansion = DataFrame(CSVFiles.load(joinpath(output_dir, "thermal_expansion.csv")))

function normalize_data!(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    if !isnothing(norm_yrs)
        idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
        for row in axes(dat, 1)
            foreach(col -> dat[row, col] -= mean(dat[row, idx_norm]), axes(dat, 2))
        end
    end
    return dat
end

normalize_data!(temperature, 1850:1900)
normalize_data!(gmslr, [2000])
normalize_data!(antarctic, [2000])
normalize_data!(greenland, [2000])

function compute_norm_quantiles(dat, norm_yrs=nothing)
    # normalize to relevant period  defined by norm_yrs
    if !isnothing(norm_yrs)
        idx_norm = findall((!isnothing).(indexin(names(dat), string.(norm_yrs))))
        for row in axes(dat, 1)
            foreach(col -> dat[row, col] -= mean(dat[row, idx_norm]), axes(dat, 2))
        end
    end
    # compute median and 95% prediction interval
    quantiles = mapcols(col -> quantile(col, [0.025, 0.5, 0.975]), dat)
    return quantiles
end

emissions_q = compute_norm_quantiles(emissions)

## load CO2 emissions for  CMIP6 scenarios
cmip_df =  XLSX.readtable("data/cmip6_co2.xlsx", "data") |> DataFrame
# select only rows and columns with scenario names and data   
select!(cmip_df, Not([:Model, :Region, :Variable, :Unit, :Notes]))
cmip_df = cmip_df[1:7, :]
# reformat scenario names to SSPn-x.x
cmip_df[!, :Scenario] = replace.(cmip_df[!, :Scenario], r"(\d)(\d)" =>   s"\1.\2")
cmip_df[!, :Scenario] = [split(cmip_df[i, :Scenario], " ")[1] for i in 1:nrow(cmip_df)]
# sort by scenarios
sort!(cmip_df, :Scenario)
# convert emissions  from MtCO2/yr to GtCO2/yr
cmip_df[!, Not(:Scenario)] = cmip_df[!, Not(:Scenario)] ./ 1000

# get empirical CDF of emissions in 2100
emis_cdf = ecdf(emissions[!, :"2100"])



fig = Figure(resolution=(1200, 800), fontsize=20, figure_padding=20)

axemissions = Axis(fig[1, 1], xlabel="Year", ylabel="CO₂ Emissions (Gt CO₂/yr)")
leg_med = Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector(emissions_q[2, :]), color=:black)
leg_ci = Makie.band!(axemissions, parse.(Int64, names(emissions_q)), Vector(emissions_q[1, :]), Vector(emissions_q[3, :]), color=(:gray, 0.2))
leg_ext = Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector((col -> maximum(col)).(eachcol(emissions))), color=:gray)
Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector((col -> minimum(col)).(eachcol(emissions))), color=:gray)
cmip_legend = Vector{LineElement}(undef, nrow(cmip_df)) # initialize storage for legend
cmip_yrs = parse.(Float64, string.(names(cmip_df[!, Not(:Scenario)])))
cmip_colors = cgrad(:Dark2_7, 7, categorical=true)
for i = 1:nrow(cmip_df)
    Makie.lines!(axemissions, cmip_yrs, Vector(cmip_df[i, Not(:Scenario)]), color=cmip_colors[i], linestyle=:dash)
    cmip_legend[i] = LineElement(color = cmip_colors[i], linestyle = :dash, linewidth=4)
end
Makie.xlims!(axemissions, 2000, 2300)
Label(fig[1, 1, TopLeft()], "a", fontsize=22, font=:bold, padding = (0, 50, 20, 0), halign=:right)
Legend(fig[2, 1], [leg_med, leg_ext, leg_ci], ["Median", "Extrema", "95% Projection Interval"], titlevisible=false, orientation=:horizontal, framevisible=true, tellheight=true)

axcdf = Axis(fig[3, 1], ylabel="CDF", xlabel="CO₂ Emissions in 2100 (Gt CO₂/yr)")
emis_range = 0:maximum(emissions[!, :"2100"])
Makie.lines!(axcdf, emis_range, emis_cdf.(emis_range), color=:black)
for i = 1:nrow(cmip_df)
    Makie.vlines!(axcdf, cmip_df[i, end], color=cmip_colors[i], linestyle=:dash)
end
Label(fig[3, 1, TopLeft()], "b", fontsize=22, font=:bold, padding = (0, 50, 20, 0), halign=:right)

leg = Legend(fig[4, 1], cmip_legend, cmip_df[!, :Scenario], "Emissions Scenario", orientation=:horizontal, tellwidth=false, tellheight=true, framevisible=false, nbanks=2)

axscatter = Axis(fig[1:2,2:3], xlabel="Global Temperature Anomaly (°C)", ylabel="Global Sea Level Anomaly (m)")
colors = cgrad(:vik100, [0.1, 0.5, 0.9])

idx2100 = findfirst(names(gmslr) .== "2100")
idx2020 = findfirst(names(gmslr) .== "2020")
cum_emissions = [sum(emissions[i, idx2020:idx2100]) for i in 1:nrow(gmslr)]
idx_low_emissions = findall(cum_emissions .< 1500)

plt =  Makie.scatter!(axscatter, temperature[:,idx2100  ], gmslr[:, idx2100], color=cum_emissions, colormap=colors, markersize=6)
cbscatter = Colorbar(fig[1:2, 4], plt, label="Cumulative CO₂ Emissions (Gt CO₂)")
Label(fig[1, 2, TopLeft()], "c", fontsize=22, font=:bold, padding = (0, 50, 20, 0), halign=:right)

axtrace = Axis(fig[3:4,2], xlabel="Year", ylabel="Global Sea Level Anomaly (m)")
for i in idx_low_emissions
    Makie.lines!(axtrace, 2020:2100, Vector(gmslr[i, idx2020:idx2100]), color=parameters[i, :t_peak], colorrange=(2030, 2050), colormap=cgrad(:Reds))
end
xlims!(axtrace, (2020, 2100))
axright = Axis(fig[3:4, 3], limits=((0, nothing), (0, 1.8)))
linkyaxes!(axtrace, axright)
Makie.density!(axright, gmslr[:, idx2100],direction=:y)
hidedecorations!(axright)
hidespines!(axright)

cbtrace = Colorbar(fig[3:4, 4], limits=(2030, 2050), colormap=:Reds, label="Year In Which Emissions Peak")
Label(fig[3, 2, TopLeft()], "d", fontsize=22, font=:bold, padding = (0, 50, 20, 0), halign=:right)


colsize!(fig.layout, 2, Auto(0.5))
colsize!(fig.layout, 3, Auto(0.1))
colgap!(fig.layout, 2, Relative(0.0))
colgap!(fig.layout, 1, Relative(0.1))

rowgap!(fig.layout, 2, Relative(0.01))
rowgap!(fig.layout, 3, Relative(0.1))

rowsize!(fig.layout, 1, Relative(0.55))
rowsize!(fig.layout, 2, Relative(0.05))
rowsize!(fig.layout, 3, Relative(0.35))
rowsize!(fig.layout, 4, Relative(0.1))

yspace = maximum(tight_yticklabel_spacing!, [axscatter, axtrace])
axscatter.yticklabelspace = yspace
axtrace.yticklabelspace = yspace


CairoMakie.save("figures/slr_temps.png", fig)





