#= this script plots scatterplots (sea level rise vs. temperature) for both
Greenland and Antarctic melting in three years: 2100, 2200, and 2300. =#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures
include("../src/functions.jl")

run = "default"
years = ["2100", "2200", "2300"]

# get the relevant results
parameters          = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "parameters.csv")))
temperature         = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "temperature.csv")))
antarctic           = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "antarctic.csv")))
greenland           = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "greenland.csv")))

# establish indices for different buckets of peaking times
early  = findall(parameters.t_peak .< 2050)           # before 2050
middle = findall(2050 .<= parameters.t_peak .< 2100)  # at least 2050 and before 2100
late   = findall(parameters.t_peak .>= 2100)          # at least 2100

# intialize relevant variables
buckets = [middle, early, late]
bucket_names = ["Middle Peaking", "Early Peaking", "Late Peaking"]
colors = [:lightblue, :seagreen4, :plum4]

# ---------------------------------- Antarctic Ice Sheet Melt vs. Temperature ---------------------------------- #

all_plots1 = []
for year in years # loop through years
    # initialize scatterplot
    p1 = scatter(title="$year", xlabel="Global Mean Temperature Anomaly (K)", ylabel="Global Mean Sea Level Contribution (m)")
    for (i,bucket) in enumerate(buckets) # loop through peaking groups
        bucket_name = bucket_names[i]
        scatter!(p1, temperature[bucket,"$year"], antarctic[bucket,"$year"], label="$bucket_name", color=colors[i],
                markerstrokewidth=0, markersize=3)
    end
    push!(all_plots1, p1) # add current year's plot to vector
end
p2 = plot(all_plots1..., plot_title="Antarctic Ice Sheet Melt vs. Temperature", layout=(1,length(years)),
        size=(600*length(years),500), ann=((0,1.1), :auto), top_margin=5mm, bottom_margin=10mm, left_margin=10mm)
#display(p2)

# ---------------------------------- Greenland Ice Sheet Melt vs. Temperature ---------------------------------- #

all_plots2 = []
for year in years # loop through years
    # initialize scatterplot
    p3 = scatter(title="$year", xlabel="Global Mean Temperature Anomaly (K)", ylabel="Global Mean Sea Level Contribution (m)")
    for (i,bucket) in enumerate(buckets) # loop through peaking groups
        bucket_name = bucket_names[i]
        scatter!(p3, temperature[bucket,"$year"], greenland[bucket,"$year"], label="$bucket_name", color=colors[i],
                markerstrokewidth=0, markersize=3)
    end
    push!(all_plots2, p3) # add current year's plot to vector
end
annotations = [((0,1.1), "(d)") ((0,1.1), "(e)") ((0,1.1), "(f)")]
p4 = plot(all_plots2..., plot_title="Greenland Ice Sheet Melt vs. Temperature", layout=(1,length(years)),
        size=(600*length(years),500), ann=annotations, top_margin=5mm, bottom_margin=10mm, left_margin=10mm)
#display(p4)

# ---------------------------------- Combine sets of plots for both GIS and AIS ---------------------------------- #

all_plots = [p2, p4]
p5 = plot(all_plots..., layout=(length(all_plots),1), size=(1500,1200))
display(p5)
savefig(p5, "/Users/ced227/Desktop/plots/icesheet_scatterplots.pdf")