# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSV of Shapley indices
using XLSX
using DataFrames # data structure for indices
using Makie # plotting library
using CairoMakie
using Measures # adjust margins with explicit measures
using StatsBase # get mean

# load the results
# Assume that we call this script from the project folder
output_dir = "output/default"
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
temperature = DataFrame(CSVFiles.load(joinpath(output_dir, "temperature.csv")))
gmslr = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
antarctic = DataFrame(CSVFiles.load(joinpath(output_dir, "antarctic.csv")))
gsic = DataFrame(CSVFiles.load(joinpath(output_dir, "gsic.csv")))
greenland = DataFrame(CSVFiles.load(joinpath(output_dir, "greenland.csv")))
lw_storage = DataFrame(CSVFiles.load(joinpath(output_dir, "lw_storage.csv")))
thermal_expansion = DataFrame(CSVFiles.load(joinpath(output_dir, "thermal_expansion.csv")))

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
temperature_q = compute_norm_quantiles(temperature, 1850:1900)
gmsl_q = compute_norm_quantiles(gmslr, [2000])
antarctic_q = compute_norm_quantiles(antarctic, [2000])
greenland_q = compute_norm_quantiles(greenland, [2000])
other_q = compute_norm_quantiles(gsic .+ lw_storage .+ thermal_expansion, [2000])

## make plot
fig = Figure(resolution=(1000, 690), fontsize=16, figure_padding=20)

gemissions = fig[1:3, 1] = GridLayout()
axemissions = Axis(gemissions[1:2, 1], xlabel="Year", ylabel="(Gt CO₂/yr)", title="Annual CO₂ Emissions", titlealign=:left, titlesize=18)
Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector(emissions_q[2, :]), color=:black)
Makie.band!(axemissions, parse.(Int64, names(emissions_q)), Vector(emissions_q[1, :]), Vector(emissions_q[3, :]), color=(:gray, 0.2))
Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector((col -> maximum(col)).(eachcol(emissions))), color=:gray)
Makie.lines!(axemissions, parse.(Int64, names(emissions_q)), Vector((col -> minimum(col)).(eachcol(emissions))), color=:gray)
cmip_legend = Vector{LineElement}(undef, nrow(cmip_df)) # initialize storage for legend
cmip_yrs = parse.(Float64, string.(names(cmip_df[!, Not(:Scenario)])))
cmip_colors = cgrad(:Dark2_7, 7, categorical=true)
for i = 1:nrow(cmip_df)
    Makie.lines!(axemissions, cmip_yrs, Vector(cmip_df[i, Not(:Scenario)]), color=cmip_colors[i], linestyle=:dash)
    cmip_legend[i] = LineElement(color = cmip_colors[i], linestyle = :dash, linewidth=4)
end
Makie.xlims!(axemissions, 2000, 2300)
leg = Legend(gemissions[3,1], cmip_legend, cmip_df[!, :Scenario], "Emissions Scenario", orientation=:horizontal, tellwidth=false, tellheight=true, framevisible=false, nbanks=4)
Label(fig[1, 1, TopLeft()], "a", fontsize=22, font=:bold, padding = (0, 35, 20, 0), halign=:right)

axtemps = Axis(fig, xlabel="Year", ylabel="Anomaly from\n1880-1900 Mean (°C)", title="Global Mean Temperature", titlealign=:left, titlesize=18, xminorticks=IntervalsBetween(4))
axgmsl = Axis(fig, xlabel="Year", ylabel="Anomaly from 2000 (m)", title="Global Mean Sea Level", titlealign=:left, titlesize=18, xminorticks=IntervalsBetween(4))
fig[1, 2] = axtemps
fig[2:3, 2] = axgmsl
linkxaxes!(axtemps, axgmsl)
# plot temperatures
Makie.lines!(axtemps, parse.(Int64, names(temperature_q)), Vector(temperature_q[2, :]), color=:black)
Makie.band!(axtemps, parse.(Int64, names(temperature_q)), Vector(temperature_q[1, :]), Vector(temperature_q[3, :]), color=(:gray, 0.2))
leg_med = Makie.lines!(axgmsl, parse.(Int64, names(gmsl_q)), Vector(gmsl_q[2, :]), color=:black)
leg_ci = Makie.band!(axgmsl, parse.(Int64, names(gmsl_q)), Vector(gmsl_q[1, :]), Vector(gmsl_q[3, :]), color=(:gray, 0.2))
Makie.xlims!(axgmsl, 2000, 2300)
hidexdecorations!(axtemps, ticks=false, grid=false, minorgrid=false)
axislegend(axgmsl, [leg_med, leg_ci], ["Median", "95% Projection Interval"], titlevisible=false, position=:lt, orientation=:vertical)
Label(fig[1, 2, TopLeft()], "b", fontsize=22, font=:bold, padding = (0, 35, 20, 0), halign=:right)
Label(fig[2, 2, TopLeft()], "c", fontsize=22, font=:bold, padding = (0, 35, 20, 0), halign=:right)

axant = Axis(fig, xlabel="Year", ylabel="(m SLE)", title="Antarctic Ice Sheet Melting", titlealign=:left, titlesize=18, xminorticks=IntervalsBetween(4))
axgreen = Axis(fig, xlabel="Year", ylabel="(m SLE)", title="Greenland Ice Sheet Melting", titlealign=:left, titlesize=18, xminorticks=IntervalsBetween(4))
axother = Axis(fig, xlabel="Year", ylabel="(m SLE)", title="Sea Level from Other Sources", titlealign=:left, titlesize=18, xminorticks=IntervalsBetween(4))
linkxaxes!(axant, axgreen, axother)
Makie.lines!(axant, parse.(Int64, names(antarctic_q)), Vector(antarctic_q[2, :]), color=:black)
Makie.band!(axant, parse.(Int64, names(antarctic_q)), Vector(antarctic_q[1, :]), Vector(antarctic_q[3, :]), color=(:gray, 0.2))
Makie.lines!(axgreen, parse.(Int64, names(greenland_q)), Vector(greenland_q[2, :]), color=:black)
Makie.band!(axgreen, parse.(Int64, names(greenland_q)), Vector(greenland_q[1, :]), Vector(greenland_q[3, :]), color=(:gray, 0.2))
Makie.lines!(axother, parse.(Int64, names(other_q)), Vector(other_q[2, :]), color=:black)
Makie.band!(axother, parse.(Int64, names(other_q)), Vector(other_q[1, :]), Vector(other_q[3, :]), color=(:gray, 0.2))
Makie.xlims!(axother, 2000, 2300)
yspace = maximum(tight_yticklabel_spacing!, [axant, axgreen, axother])
axant.yticklabelspace = yspace
axgreen.yticklabelspace = yspace
axother.yticklabelspace = yspace
hidexdecorations!(axant, ticks=false, grid=false, minorgrid=false)
hidexdecorations!(axgreen, ticks=false, grid=false, minorgrid=false)
fig[1, 3] = axant
fig[2, 3] = axgreen
fig[3, 3] = axother
Label(fig[1, 3, TopLeft()], "d", fontsize=22, font=:bold, padding = (0, 35, 20, 0), halign=:right)
Label(fig[2, 3, TopLeft()], "e", fontsize=22, font=:bold, padding = (0, 35, 20, 0), halign=:right)
Label(fig[3, 3, TopLeft()], "f", fontsize=22, font=:bold, padding = (0, 35, 20, 0), halign=:right)

CairoMakie.save("figures/ensemble_projections.png", fig)
