#= CHANGE this script divides the results of model_ensemble.jl into three buckets based on peaking time: before 2050, 2050-2100, and after 2100.
Then plots boxplots for relevant outputs for each of the three buckets for three years: 2100, 2200, and 2300. 
Plot 3: shows the corresponding GMSLR contributions for ONE run for each of low,med,high emissions (doesn't show uncertainty) =#

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
quantiles   = mapslices(x -> quantile(x, [0.05, 0.475, 0.525, 0.95]), Matrix(gmslr), dims=1) # 5%, 47.5%, 52.5%, and 95% CI
low_gmslr   = findall(gmslr[:,end] .< quantiles[1,end])                                          # lowest 5% of GMSLR outcomes
med_gmslr   = findall((gmslr[:,end] .> quantiles[2,end]) .&& (gmslr[:,end] .< quantiles[3,end])) # middle 5% of GMSLR outcomes
high_gmslr  = findall(gmslr[:,end] .> quantiles[4,end])                                          # highest 5% of GMSLR outcomes

# define the index for the scenarios that will be plotted
Random.seed!(1)
s1 = rand(intersect(early_tpeak,  low_gmslr))
s2 = rand(intersect(middle_tpeak, low_gmslr))
s3 = rand(intersect(late_tpeak,   low_gmslr))
s4 = rand(intersect(early_tpeak,  med_gmslr))
s5 = rand(intersect(middle_tpeak, med_gmslr))
s6 = rand(intersect(late_tpeak,   med_gmslr))
s7 = rand(intersect(early_tpeak,  high_gmslr))
s8 = rand(intersect(middle_tpeak, high_gmslr))
s9 = rand(intersect(late_tpeak,   high_gmslr))

# -------------------------------------------------------------------------------------------------------- #
# --------------- Plot: Comparison of GMSLR Contributors for Different Scenarios ------------------------- #
# -------------------------------------------------------------------------------------------------------- #

all_plots = []
scenarios = [s1, s2, s3, s4, s5, s6, s7, s8, s9]
# define titles for each plot
titles = ["Early Peaking, Low GMSLR",    "Middle Peaking, Low GMSLR",    "Late Peaking, Low GMSLR",
          "Early Peaking, Medium GMSLR", "Middle Peaking, Medium GMSLR", "Late Peaking, Medium GMSLR",
          "Early Peaking, High GMSLR",   "Middle Peaking, High GMSLR",   "Late Peaking, High GMSLR"]
# define labels for plotting
labels = ["Antarctic" "GSIC" "Greenland" "LW Storage" "Thermal Expansion"]
years = parse.(Int64, names(gmslr)) # vector of years
xticks = first(years):50:last(years)
ylims = [(-0.9,1.2), (-0.4,4), (-0.5,8)]

for (i,scenario) in enumerate(scenarios) # loop through scenarios
    j = trunc(Int64, (i-1)/3)+1 # index for yticks
    # initialize the plot
    p = plot(title=titles[i], ylim=ylims[j], xlabel="Year", ylabel="Global Mean Sea Level Anomaly (m)")
    # create a matrix of sea level rise contributors
    slr_contributors = [collect(antarctic[scenario,:])';         # Antarctic Ice Sheet
                        collect(gsic[scenario,:])';              # glaciers and small ice sheets
                        collect(greenland[scenario,:])';         # Greenland Ice Sheet
                        collect(lw_storage[scenario,:])';        # land water storage
                        collect(thermal_expansion[scenario,:])'] # thermal expansion
    areaplot!(years, slr_contributors', label=labels, color_palette=:darkrainbow, alpha=1, fillalpha=0.5)
    plot!(years, collect(gmslr[scenario,:]), label="GMSLR", color=:black, linewidth=4, linestyle=:dash)
    # add current plot to vector
    push!(all_plots, p)
end

# combine all plots and save
p2 = plot(all_plots..., size=(1800,1200), xticks=xticks, bottom_margin=7mm, left_margin=8mm, legend=:topleft)
display(p2)
#savefig(p2, "/Users/ced227/Desktop/plots/scenario_combos_stacked_area.png")

# -------------------------------------------------------------------------------------------------------- #
# --------------------- Supplementary: Can be useful for benchmarking ------------------------------------ #
# -------------------------------------------------------------------------------------------------------- #

#= subset df to just include emissions parameters
emission_params = parameters[:,[:gamma_g, :t_peak, :gamma_d]]

# establish relevant quantiles to divide emissions trajectories
q_middle = mapslices(x -> quantile(x, [0.45, 0.55]), Matrix(emission_params), dims=1) # middle 10% of samples
q_80     = mapslices(x -> quantile(x, [0.1, 0.9]), Matrix(emission_params), dims=1) # 80% credible interval
q_90     = mapslices(x -> quantile(x, [0.05, 0.95]), Matrix(emission_params), dims=1) # 90% credible interval

# isolate column for each emissions parameter
growth = parameters[:,:gamma_g]
peak = parameters[:,:t_peak]
decline = parameters[:,:gamma_d]

# initialize storage for indices translating to low, medium, and high emissions samples
lower_idx = []
medium_idx = []
upper_idx = []

for i in 1:num_samples
    # if we have early peaking and rapid decarbonization
    if peak[i] <= q_90[1,2] && decline[i] >= q_90[2,3]
        push!(lower_idx, i)
    # if we have medium peaking and medium decarbonization
    elseif q_middle[1,2] <= peak[i] <= q_middle[2,2] && q_middle[1,3] <= decline[i] <= q_middle[2,3]
        push!(medium_idx, i)
    # if we have late peaking and slow decarbonization
    elseif peak[i] >= q_90[2,2] && decline[i] <= q_90[1,3]
        push!(upper_idx, i)
    end
end
println(length(lower_idx)) # shows how many samples meet the criteria
println(length(medium_idx))
println(length(upper_idx))

# df of values for each emission parameter meeting criteria (growth, peak, decline) (for reference)
lower_samples  = emission_params[lower_idx,:]
medium_samples = emission_params[medium_idx,:]
upper_samples  = emission_params[upper_idx,:]
=#