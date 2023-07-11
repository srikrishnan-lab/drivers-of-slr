#= this script divides the results of model_ensemble.jl into scenario combinations based on three buckets of GMSLR outcomes and three buckets of peaking time.
Then plots stacked area plots for ONE sample (doesn't show uncertainty) for each combination of peaking time & GMSLR (9 scenario combinations total). =#

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

# establish indices for different buckets of GMSLR
quantiles  = mapslices(x -> quantile(x, [0.05, 0.475, 0.525, 0.95]), Matrix(gmslr), dims=1) # 5%, 47.5%, 52.5%, and 95% CI
low_gmslr  = findall(gmslr[:,end] .< quantiles[1,end])                                          # lowest 5% of GMSLR outcomes
med_gmslr  = findall((gmslr[:,end] .> quantiles[2,end]) .&& (gmslr[:,end] .< quantiles[3,end])) # middle 5% of GMSLR outcomes
high_gmslr = findall(gmslr[:,end] .> quantiles[4,end])                                          # highest 5% of GMSLR outcomes

# establish indices for different buckets of peaking times
early_tpeak  = findall(parameters.t_peak .< 2050)           # before 2050
middle_tpeak = findall(2050 .<= parameters.t_peak .< 2100)  # at least 2050 and before 2100
late_tpeak   = findall(parameters.t_peak .>= 2100)          # at least 2100

# establish the indices for samples that fit the criteria for each scenario combination
s1 = intersect(low_gmslr,  early_tpeak)
s2 = intersect(low_gmslr,  middle_tpeak)
s3 = intersect(low_gmslr,  late_tpeak)
s4 = intersect(med_gmslr,  early_tpeak)
s5 = intersect(med_gmslr,  middle_tpeak)
s6 = intersect(med_gmslr,  late_tpeak)
s7 = intersect(high_gmslr, early_tpeak)
s8 = intersect(high_gmslr, middle_tpeak)
s9 = intersect(high_gmslr, late_tpeak)

# create a dataframe that shows the number of samples that fit each scenario's criteria
scenarios = [s1, s2, s3, s4, s5, s6, s7, s8, s9]
scenario_names = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9"]
sea_level = ["Low", "Low", "Low", "Medium", "Medium", "Medium", "High", "High", "High"]
peaking = ["Early", "Middle", "Late", "Early", "Middle", "Late", "Early", "Middle", "Late"]
df = DataFrame("Scenario Combination" => scenario_names, "GMSLR Group" => sea_level, "Peaking Group" => peaking, "Number of Samples" => length.(scenarios))
#XLSX.writetable("/Users/ced227/Desktop/plots/scenario_combos.xlsx", "scenarios" => df) # save df to Excel

# randomly select one index from vector of samples that fit criteria for plotting
Random.seed!(1)
s1 = rand(s1)
s2 = rand(s2)
s3 = rand(s3)
s4 = rand(s4)
s5 = rand(s5)
s6 = rand(s6)
s7 = rand(s7)
s8 = rand(s8)
s9 = rand(s9)

# -------------------------------------------------------------------------------------------------------- #
# --------------- Plot: Comparison of GMSLR Contributors for Different Scenarios ------------------------- #
# -------------------------------------------------------------------------------------------------------- #

all_plots = []
scenarios = [s1, s2, s3, s4, s5, s6, s7, s8, s9]
years = parse.(Int64, names(gmslr)) # vector of years
xticks = first(years):50:last(years)
ylims = [(-0.9,1.2), (-0.4,4), (-0.5,8)]
labels = ["Antarctic" "GSIC" "Greenland" "LW Storage" "Thermal Expansion"] # labels for plot legend
# define titles for each plot
titles = ["S1: Low GMSLR, Early Peaking",    "S2: Low GMSLR, Middle Peaking",    "S3: Low GMSLR, Late Peaking",
          "S4: Medium GMSLR, Early Peaking", "S5: Medium GMSLR, Middle Peaking", "S6: Medium GMSLR, Late Peaking",
          "S7: High GMSLR, Early Peaking",   "S8: High GMSLR, Middle Peaking",   "S9: High GMSLR, Late Peaking"]

for (i,scenario) in enumerate(scenarios) # loop through scenarios
    j = trunc(Int64, (i-1)/3)+1 # index for yticks
    # initialize the plot
    p1 = plot(title=titles[i], ylim=ylims[j], xlabel="Year", ylabel="Global Mean Sea Level Anomaly (m)")
    # create a matrix of sea level rise contributors
    slr_contributors = [collect(antarctic[scenario,:])';         # Antarctic Ice Sheet
                        collect(gsic[scenario,:])';              # glaciers and small ice sheets
                        collect(greenland[scenario,:])';         # Greenland Ice Sheet
                        collect(lw_storage[scenario,:])';        # land water storage
                        collect(thermal_expansion[scenario,:])'] # thermal expansion
    areaplot!(years, slr_contributors', label=labels, color_palette=:darkrainbow, alpha=1, fillalpha=0.5)
    plot!(years, collect(gmslr[scenario,:]), label="GMSLR", color=:black, linewidth=4, linestyle=:dash)
    # add current plot to vector
    push!(all_plots, p1)
end

# combine all plots and save
p2 = plot(all_plots..., size=(1800,1200), xticks=xticks, bottom_margin=7mm, left_margin=8mm, legend=:topleft, ann=((0,1.06), :auto))
display(p2)
#savefig(p2, "/Users/ced227/Desktop/plots/scenario_combos_stacked_area.pdf")

# -------------------------------------------------------------------------------------------------------- #
# ---------------- Analysis of Parameter/Outcome Percentiles for each Scenario --------------------------- #
# -------------------------------------------------------------------------------------------------------- #
let
    cum_emissions = sum.(eachrow(emissions)) # cumulative emissions in 2300 for all samples
    percentile_emissions = zeros(length(scenarios)) # percentiles for cumulative emissions
    percentile_temps = zeros(length(scenarios)) # percentiles for temperature anomaly in 2300

    # create df and add necessary columns
    param_names = names(parameters)
    all_percentiles = DataFrame([scenario => zeros(length(param_names)) for scenario in scenario_names])
    insertcols!(all_percentiles, 1, :parameter => param_names)

    for (i, scenario) in enumerate(scenarios) # loop through scenario indices
        idx = scenario # establish index of sample we want to consider

        for (j, param) in enumerate(param_names) # loop through parameters
            # find the value and percentile for the sample
            value = parameters[idx, param]
            percentile = 100 * mean(parameters[:, param] .< value)
            # add the percentile to df
            all_percentiles[j,i+1] = percentile
        end

        # include percentile for cumulative emissions
        percentile_emissions[i] = 100 * mean(cum_emissions .< cum_emissions[idx])
        # include percentile for temperature in 2300
        percentile_temps[i] = 100 * mean(temperature[:,end] .< temperature[:,end][idx])

        if i == length(scenarios) # if we're on the last run
            # insert emissions percentiles into first row of all_percentiles df 
            insert!.(eachcol(all_percentiles), 1, ["cumulative_emissions", percentile_emissions...])
            # insert temperature percentiles into first row of all_percentiles df 
            insert!.(eachcol(all_percentiles), 1, ["temperature", percentile_temps...])
        end
    end
    show(all_percentiles, allrows=true)
end
# get cumulative emissions or temperature for a certain sample
#cum_emissions[s1]
#temperature[:,end][s1]