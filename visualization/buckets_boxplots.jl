#= this script divides the results of model_ensemble.jl into three buckets based on peaking time: before 2050, 2050-2100, and after 2100.
Then plots boxplots for each of the three buckets for three years: 2100, 2200, and 2300. =#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures, StatsPlots
include("../src/functions.jl")

output_dir = "default"

# load the results
parameters          = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "parameters.csv")))
emissions           = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "emissions.csv")))
radiative_forcing   = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "radiative_forcing.csv")))
temperature         = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "temperature.csv")))
gmslr               = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "gmslr.csv")))
antarctic           = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "antarctic.csv")))
gsic                = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "gsic.csv")))
greenland           = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "greenland.csv")))
lw_storage          = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "lw_storage.csv")))
thermal_expansion   = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "thermal_expansion.csv")))
ocean_heat          = DataFrame(load(joinpath(@__DIR__, "..", "results", "$output_dir", "ocean_heat.csv")))

# establish indices for different buckets of peaking times
early_indices   = findall(parameters.t_peak .<= 2050)           # at or before 2050
middle_indices  = findall(2050 .< parameters.t_peak .< 2100)    # between 2050 and 2100
late_indices    = findall(parameters.t_peak .>= 2100)           # at or after 2100

let
    all_years = []
    years = ["2100", "2200", "2300"]
    buckets = [early_indices, middle_indices, late_indices]
    #loop through years
    for (i, year) in enumerate(years)
        all_buckets = []
        # loop through buckets of early, middle, and late peaking
        for (j, bucket) in enumerate(buckets)
            # boxplots for rf, temp, and GMSLR
            xticks = (1:3, ["Radiative Forcing" "Temperature" "GMSLR"])
            ylims = [(0,9), (0,11), (0,13)]
            p1 = boxplot(radiative_forcing[bucket,"$year"], color=:purple, xticks=xticks, ylims=ylims[i])
            boxplot!(temperature[bucket,"$year"], color=:grey)
            boxplot!(gmslr[bucket,"$year"], color=:red)
            # boxplots for contributors to GMSLR
            xticks = (1:5, ["Antarctic" "GSIC" "Greenland" "LW Storage" "Thermal Expansion"])
            ylims = [(-0.1, 0.9), (-0.4, 3), (-0.6, 5)]
            p2 = boxplot(antarctic[bucket,"$year"], xticks=xticks, ylims=ylims[i])
            boxplot!(gsic[bucket,"$year"])
            boxplot!(greenland[bucket,"$year"])
            boxplot!(lw_storage[bucket,"$year"])
            boxplot!(thermal_expansion[bucket,"$year"])
            # combined plots of p1 and p2
            titles = ["Early Peaking: $year" "Middle Peaking: $year" "Late Peaking: $year"]
            p3 = plot(p1, p2, layout=(1,2), size=(700,400), xrotation=45, legend=:false, plot_title=titles[j])
            push!(all_buckets, p3)
        end 
        # combined plots for early, middle, and late peaking
        p4 = plot(all_buckets..., layout=(1,3), size=(2100,400), legend=:false)
        push!(all_years, p4)
    end
    # combined plots for all years
    p5 = plot(all_years..., layout=(3,1), size=(2100,1200))
    display(p5)
    #savefig(p5, "/Users/ced227/Desktop/plots/boxplots.png")
end