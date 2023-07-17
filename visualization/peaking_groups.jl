# this script creates three plots, one for each peaking group, to visualize groups for defense presentation
# each plot includes three selected trajectories for each group and a shaded region for the lower/upper bound of peaking time

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures
include("../src/functions.jl")

run_name = "default"

# get the relevant results
parameters = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run_name", "parameters.csv")))
emissions  = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run_name", "emissions.csv")))

# establish indices for different buckets of peaking times
early  = findall(parameters.t_peak .< 2050)           # before 2050
middle = findall(2050 .<= parameters.t_peak .< 2100)  # at least 2050 and before 2100
late   = findall(parameters.t_peak .>= 2100)          # at least 2100

# initialize historical observations
historical_data = historical_emissions() # years 1850-2021

# initialize values
start_year = 2022
end_year = 2300
peaking_groups = [early, middle, late]
group_names = ["Early", "Middle", "Late"]

# get 10%, median, and 90% quantiles for growth and decline (same for all groups)
growth_quantiles = quantile(parameters[:, "gamma_g"], (0.1, 0.5, 0.9))
decline_quantiles = quantile(parameters[:, "gamma_d"], (0.1, 0.5, 0.9))

# create arrays with three levels of growth/decline rates
growth = [growth_quantiles[1], growth_quantiles[2], growth_quantiles[3]] # slowest to fastest growth
decline = [decline_quantiles[3], decline_quantiles[2], decline_quantiles[1]] # fastest to slowest decline

# loop through peaking groups
for (i, group) in enumerate(peaking_groups)

    group_name = group_names[i] # establish name for current group
    p = plot(title="$group_name Peaking Trajectories", xlabel="Year", ylabel="CO₂ Emissions (GtCO₂/yr)", 
            xticks=collect(2000:50:end_year), ylims=(0,130), palette=:seaborn_dark)

    # get 10%, median, and 90% quantiles for current group's peaking time
    peaking_quantiles = quantile(parameters[group, "t_peak"], (0.1, 0.5, 0.9))

    # create three representative emissions curves for current group (incremental growth, peak, decline for each curve)
    t1, q1 = emissions_curve(historical_data, γ_g=growth[1], t_peak=peaking_quantiles[1], γ_d=decline[1], start_year=start_year, end_year=end_year)
    t2, q2 = emissions_curve(historical_data, γ_g=growth[2], t_peak=peaking_quantiles[2], γ_d=decline[2], start_year=start_year, end_year=end_year)
    t3, q3 = emissions_curve(historical_data, γ_g=growth[3], t_peak=peaking_quantiles[3], γ_d=decline[3], start_year=start_year, end_year=end_year)

    plot!(p, t1, q1, linewidth=4, label="10th Percentile")
    plot!(p, t2, q2, linewidth=4, label="Median")
    plot!(p, t3, q3, linewidth=4, label="90th Percentile")

    # add historical points from 2000-2021
    indices = findall((in)(2000:last(historical_data.year)), first(historical_data.year):last(historical_data.year))
    scatter!(historical_data.year[indices], historical_data.rcp_co2_emissions[indices], label="Historical Data", markersize=3, color=:black)

    # shade region for peaking group
    if i == 1 # if early peaking group
        vspan!([2030, 2049], linecolor=false, fillcolor=:grey, fillalpha=0.5, label="Peaking Group Bounds")
    elseif i == 2 # if middle peaking group
        vspan!([2050, 2099], linecolor=false, fillcolor=:grey, fillalpha=0.5, label="Peaking Group Bounds")
    elseif i == 3 # if late peaking group
        vspan!([2100, 2199], linecolor=false, fillcolor=:grey, fillalpha=0.5, label="Peaking Group Bounds")
    end

    # change font size and display/save figure
    plot!(p, labelfontsize=13, tickfontsize=12, titlefontsize=18, legend=:false)
    display(p)
    #savefig(p, "/Users/ced227/Desktop/defense/$group_name.pdf")
end