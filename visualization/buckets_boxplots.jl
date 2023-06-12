#= this script divides the results of model_ensemble.jl into three buckets based on peaking time: before 2050, 2050-2100, and after 2100.
Then plots boxplots for relevant outputs for each of the three buckets for different years (either on centennial or decadal time scales). =#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures, StatsPlots
include("../src/functions.jl")

output_dir = "default"

# uncomment one of the two vectors below to select years
#years = ["2100", "2200", "2300"]
#years = ["2070", "2100", "2120", "2150"]

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
early_indices   = findall(parameters.t_peak .< 2050)           # before 2050
middle_indices  = findall(2050 .<= parameters.t_peak .< 2100)  # at least 2050 and before 2100
late_indices    = findall(parameters.t_peak .>= 2100)          # at least 2100

# --------------------------------------------------------------------------------- #
# -------------- Plot radiative forcing, temperature, and GMSLR ------------------- #
# --------------------------------------------------------------------------------- #
let
    all_outputs = []
    outputs = [radiative_forcing, temperature, gmslr]
    # loop through outputs of rf, temp, and GMSLR
    for (i, output) in enumerate(outputs)
        all_years = []
        # loop through years
        for year in years
            # boxplots for early, middle, and late peaking for current year
            titles = ["Radiative Forcing: $year" "Temperature: $year" "Global Mean Sea Level: $year"]
            ylabels = ["Global Radiative Forcing (W/m²)", "Global Mean Temperature Anomaly (K)", "Global Mean Sea Level Anomaly (m)"]
            xticks = (1:3, ["Early" "Middle" "Late"])
            ylims = []
            if length(years) == 3
                ylims = [(0, 15.5), (0, 11), (-0.5, 13.5)] # ylims for 2100, 2200, 2300
            elseif length(years) == 4
                ylims = [(0, 13), (0, 8.5), (-0.3, 5.5)] # ylims for 2070, 2100, 2120, 2150
            end
            p1 = boxplot(xticks=xticks, ylim=ylims[i], title=titles[i], ylabel=ylabels[i])
            boxplot!(output[early_indices,"$year"], color=:seagreen4) # early peaking
            boxplot!(output[middle_indices,"$year"], color=:lightblue) # middle peaking
            boxplot!(output[late_indices,"$year"], color=:plum4) # late peaking
            # add current year's plot to vector
            push!(all_years, p1)
        end 
        # combined plots for all years
        p2 = plot(all_years..., layout=(1,length(years)), legend=:false)
        push!(all_outputs, p2)
    end
    # combined plots for all outputs
    p3 = plot(all_outputs..., layout=(3,1), size=(1800,1200), bottom_margin=5mm, left_margin=10mm, ann=((0,1.06),:auto), dpi=600)
    display(p3)
    #savefig(p3, "/Users/ced227/Desktop/plots/rf_temp_gmslr_boxplots.png")
end

# --------------------------------------------------------------------------------- #
# ------------------------- Plot contributors to GMSLR ---------------------------- #
# --------------------------------------------------------------------------------- #
let
    all_years = []
    buckets = [early_indices, middle_indices, late_indices]
    # loop through years
    for (i, year) in enumerate(years)
        all_buckets = []
        # loop through buckets of early, middle, and late peaking
        for (j, bucket) in enumerate(buckets)
            # boxplots for contributors to GMSLR
            titles = ["Early Peaking: $year" "Middle Peaking: $year" "Late Peaking: $year"]
            xticks = (1:5, ["Antarctic" "GSIC" "Greenland" "LW Storage" "TE"])
            ylims = []
            if length(years) == 3
                ylims = [(-0.4, 1.75), (-0.9, 5.3), (-1.2, 7.2)] # ylims for 2100, 2200, 2300
            elseif length(years) == 4
                ylims = [(-0.3, 0.9), (-0.4, 1.75), (-0.5, 2.4), (-0.8, 3.5)] # ylims for 2070, 2100, 2120, 2150
            end
            p1 = boxplot(xticks=xticks, ylims=ylims[i], title=titles[j], ylabel="Global Mean Sea Level Anomaly (m)")
            boxplot!(antarctic[bucket,"$year"])         # Antarctic ice sheet
            boxplot!(gsic[bucket,"$year"])              # glaciers and small ice caps
            boxplot!(greenland[bucket,"$year"])         # Greenland ice sheet
            boxplot!(lw_storage[bucket,"$year"])        # land water storage
            boxplot!(thermal_expansion[bucket,"$year"]) # thermal expansion
            # add current bucket's plot to vector
            push!(all_buckets, p1)
        end 
        # combined plots for early, middle, and late peaking
        p2 = plot(all_buckets..., layout=(1,3), legend=:false)
        push!(all_years, p2)
    end
    # combined plots for all years
    loc = length(years)==3 ? (0,1.06) : (0,1.09) # location of panel labels
    p3 = plot(all_years..., layout=(length(years),1), size=(1800,1200), bottom_margin=5mm, left_margin=10mm, ylabelfontsize=9, ann=(loc,:auto), dpi=600)
    display(p3)
    #savefig(p3, "/Users/ced227/Desktop/plots/slr_boxplots.png")
end