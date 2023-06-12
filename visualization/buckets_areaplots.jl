#= this script divides the results of model_ensemble.jl into three buckets of peaking time and three buckets of GMSLR outcomes.
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

# establish indices for different buckets of peaking times
early_tpeak  = findall(parameters.t_peak .< 2050)           # before 2050
middle_tpeak = findall(2050 .<= parameters.t_peak .< 2100)  # at least 2050 and before 2100
late_tpeak   = findall(parameters.t_peak .>= 2100)          # at least 2100

# establish indices for different buckets of GMSLR
quantiles  = mapslices(x -> quantile(x, [0.05, 0.475, 0.525, 0.95]), Matrix(gmslr), dims=1) # 5%, 47.5%, 52.5%, and 95% CI
low_gmslr  = findall(gmslr[:,end] .< quantiles[1,end])                                          # lowest 5% of GMSLR outcomes
med_gmslr  = findall((gmslr[:,end] .> quantiles[2,end]) .&& (gmslr[:,end] .< quantiles[3,end])) # middle 5% of GMSLR outcomes
high_gmslr = findall(gmslr[:,end] .> quantiles[4,end])                                          # highest 5% of GMSLR outcomes

# establish the indices for samples that fit the criteria for each scenario combination
s1 = intersect(early_tpeak,  low_gmslr)
s2 = intersect(middle_tpeak, low_gmslr)
s3 = intersect(late_tpeak,   low_gmslr)
s4 = intersect(early_tpeak,  med_gmslr)
s5 = intersect(middle_tpeak, med_gmslr)
s6 = intersect(late_tpeak,   med_gmslr)
s7 = intersect(early_tpeak,  high_gmslr)
s8 = intersect(middle_tpeak, high_gmslr)
s9 = intersect(late_tpeak,   high_gmslr)

# create a dataframe that shows the number of samples that fit each scenario's criteria
scenarios = [s1, s2, s3, s4, s5, s6, s7, s8, s9]
peaking = ["Early", "Middle", "Late", "Early", "Middle", "Late", "Early", "Middle", "Late"]
sea_level = ["Low", "Low", "Low", "Medium", "Medium", "Medium", "High", "High", "High"]
df = DataFrame("Scenario" => collect(1:9), "Peaking Group" => peaking, "GMSLR Group" => sea_level, "Number of Samples" => length.(scenarios))
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
titles = ["Early Peaking, Low GMSLR",    "Middle Peaking, Low GMSLR",    "Late Peaking, Low GMSLR",
          "Early Peaking, Medium GMSLR", "Middle Peaking, Medium GMSLR", "Late Peaking, Medium GMSLR",
          "Early Peaking, High GMSLR",   "Middle Peaking, High GMSLR",   "Late Peaking, High GMSLR"]

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