using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using DataFrames
using CSVFiles
using Plots
using StatsPlots
using Distributions

include("functions.jl")

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

num_samples = size(parameters, 1)
i = 1 #rand(1:num_samples) # random run to look at (between 1 and num_samples)

# -------------------------------------------------------------------------------------------------------- #
# ---------------------- Plot 1: Input/Output Visualization for Individual Run --------------------------- #
# -------------------------------------------------------------------------------------------------------- #

# create variables for years and xticks for plotting
years = parse.(Int64, names(emissions))
xticks = first(years):50:last(years)

# plot of emissions: values for current run, as well as 95% CI
p1 = plot(title="CO₂ Emissions", xlabel="Year", ylabel="Total CO₂ Emissions (GtCO₂/yr)", legend=:topleft, xticks=xticks, ylim=(0,100))
plot!(years, collect(emissions[i,:]), label="Emissions", color=:blue, linewidth=2)
co2_quantiles = mapslices(x -> quantile(x, [0.025, 0.975]), Matrix(emissions), dims=1) # 95% credible interval
plot!(years, co2_quantiles[1,:], fillrange=co2_quantiles[2,:], fillalpha=0.4, alpha= 0.35, color=:blue, label="95% CI")

# plot of temperature: values for current run, as well as 95% CI
p2 = plot(title="Temperature", xlabel="Year", ylabel="Global Mean Temperature Anomaly (K)", legend=:topleft, xticks=xticks, ylim=(-1,9))
plot!(years, collect(temperature[i,:]), label="Temperature", color=:green, linewidth=2)
temp_quantiles = mapslices(x -> quantile(x, [0.025, 0.975]), Matrix(temperature), dims=1) # 95% credible interval
plot!(years, temp_quantiles[1,:], fillrange=temp_quantiles[2,:], fillalpha=0.4, alpha= 0.35, color=:green, label="95% CI")

# plot of radiative forcing: values for current run, as well as 95% CI
p3 = plot(title="Radiative Forcing", xlabel="Year", ylabel="Global Radiative Forcing (W/m²)", legend=:topleft, xticks=xticks, ylim=(-2,12))
plot!(years, collect(radiative_forcing[i,:]), label="Radiative Forcing", color=:purple, linewidth=2)
rf_quantiles = mapslices(x -> quantile(x, [0.025, 0.975]), Matrix(radiative_forcing), dims=1) # 95% credible interval
plot!(years, rf_quantiles[1,:], fillrange=rf_quantiles[2,:], fillalpha=0.4, alpha=0.35, color=:purple, label="95% CI")

# create a stacked area plot for sea level rise
p4 = plot(title="GMSLR Contributions", xlabel="Year", ylabel="Global Mean Sea Level Anomaly (m)", legend=:topleft, xticks=xticks, ylim=(-1,10))
# create a matrix of sea level rise contributors
slr_contributions = [collect(antarctic[i,:])';
                     collect(gsic[i,:])';
                     collect(greenland[i,:])';
                     collect(lw_storage[i,:])';
                     collect(thermal_expansion[i,:])'] # 5 x num_years matrix
labels=["Antarctic" "GSIC" "Greenland" "LW Storage" "Thermal Expansion"]
areaplot!(years, slr_contributions', label=labels, color_palette=:darkrainbow, alpha=1, fillalpha=0.5)
plot!(years, collect(gmslr[i,:]), label="GMSLR", color=:black, linewidth=4, linestyle=:dash)

# combine all plots and format
all_plots1 = plot(p1, p2, p3, p4, layout=4, plot_title="Run $i Results - $output_dir")
plot!(size=(1000,800), labelfontsize=10, margin=5Plots.mm)
display(all_plots1)
savefig(all_plots1, "/Users/ced227/Desktop/plots/run$i-$output_dir.png")

# -------------------------------------------------------------------------------------------------------- #
# ---- Plot 2: Representative Low, Medium, & High Emissions Trajectories with SLR Range Visualization ---- #
# -------------------------------------------------------------------------------------------------------- #

# get results for HIGH emissions and associated output
emissions_high = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "emissions.csv")))
rf_high        = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "radiative_forcing.csv")))
temp_high      = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "temperature.csv")))
gmslr_high     = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "gmslr.csv")))
antarctic_high = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "antarctic.csv")))
gsic_high      = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "gsic.csv")))
greenland_high = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "greenland.csv")))
lw_high        = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "lw_storage.csv")))
te_high        = DataFrame(load(joinpath(@__DIR__, "..", "results", "high", "thermal_expansion.csv")))
# get results for MEDIUM emissions and associated output
emissions_med  = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "emissions.csv")))
rf_med         = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "radiative_forcing.csv")))
temp_med       = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "temperature.csv")))
gmslr_med      = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "gmslr.csv")))
antarctic_med  = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "antarctic.csv")))
gsic_med       = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "gsic.csv")))
greenland_med  = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "greenland.csv")))
lw_med         = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "lw_storage.csv")))
te_med         = DataFrame(load(joinpath(@__DIR__, "..", "results", "medium", "thermal_expansion.csv")))
# get results for LOW emissions and associated output
emissions_low  = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "emissions.csv")))
rf_low         = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "radiative_forcing.csv")))
temp_low       = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "temperature.csv")))
gmslr_low      = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "gmslr.csv")))
antarctic_low  = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "antarctic.csv")))
gsic_low       = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "gsic.csv")))
greenland_low  = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "greenland.csv")))
lw_low         = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "lw_storage.csv")))
te_low         = DataFrame(load(joinpath(@__DIR__, "..", "results", "low", "thermal_expansion.csv")))

# create plots for emissions and GMSLR
plt1 = plot(title="Representative CO₂ Emissions Trajectories", xlabel="Year", ylabel="Total CO₂ Emissions (GtCO₂/yr)", legend=:topleft, xticks=xticks)
plt2 = plot(title="GMSLR Credible Intervals", xlabel="Year", ylabel="Global Mean Sea Level Anomaly (m)", legend=:topleft, xticks=xticks)

# add representative high emissions curve
plot!(plt1, years, collect(emissions_high[1,:]), label="High Emissions", color=:purple, xticks=xticks, linewidth=3)
# add representative medium emissions curve
plot!(plt1, years, collect(emissions_med[1,:]), label="Medium Emissions", color=:blue, xticks=xticks, linewidth=3)
# add representative low emissions curve
plot!(plt1, years, collect(emissions_low[1,:]), label="Low Emisions", color=:green, xticks=xticks, linewidth=3)

# add in extreme RCP scenario emissions (GtCO₂)
rcp26, rcp85 = rcp_emissions()
scatter!(plt1, rcp85[:,1], rcp85[:,2], label="RCP 8.5", color=:black, markersize=3, shape=:rect)
scatter!(plt1, rcp26[:,1], rcp26[:,2], label="RCP 2.6", color=:black, markersize=4, shape=:utriangle)

# create quantiles for low, medium, and high emissions samples for sea level rise (95% credible interval)
slr_quantiles_high = mapslices(x -> quantile(x, [0.025, 0.975]), Matrix(gmslr_high), dims=1)
slr_quantiles_med  = mapslices(x -> quantile(x, [0.025, 0.975]), Matrix(gmslr_med), dims=1)
slr_quantiles_low  = mapslices(x -> quantile(x, [0.025, 0.975]), Matrix(gmslr_low), dims=1)

# add in credible intervals for low, medium, and high sea level rise
plot!(plt2, years, slr_quantiles_high[1,:], fillrange=slr_quantiles_high[2,:], fillalpha=0.5, alpha=1, color=:purple, label="95% CI High Emissions")
plot!(plt2, years, slr_quantiles_med[1,:], fillrange=slr_quantiles_med[2,:], fillalpha=0.5, alpha=1, color=:blue, label="95% CI Medium Emissions")
plot!(plt2, years, slr_quantiles_low[1,:], fillrange=slr_quantiles_low[2,:], fillalpha=0.5, alpha=1, color=:green, label="95% CI Low Emissions")

# combine plots and display
all_plots2 = plot(plt1, plt2, layout=(1,2), size=(1100,400), margin=5Plots.mm)
display(all_plots2)
savefig(all_plots2, "/Users/ced227/Desktop/plots/low_med_high_co2_slr.png")

# -------------------------------------------------------------------------------------------------------- #
# --------------- Plot 3: Output Visualization Low, Medium, High Emissions Comparison -------------------- #
# -------------------------------------------------------------------------------------------------------- #

# plot of temperature: values for low, med, high emissions
f1 = plot(title="Temperature", xlabel="Year", ylabel="Global Mean Temperature Anomaly (K)", legend=:topleft, xticks=xticks, ylim=(-1,9))
plot!(years, collect(temp_high[1,:]), label="High Emissions", color=:purple, linewidth=2)
plot!(years, collect(temp_med[1,:]), label="Medium Emissions", color=:blue, linewidth=2)
plot!(years, collect(temp_low[1,:]), label="Low Emissions", color=:green, linewidth=2)

# plot of radiative forcing: values for low, med, high emissions
f2 = plot(title="Radiative Forcing", xlabel="Year", ylabel="Global Radiative Forcing (W/m²)", legend=:topleft, xticks=xticks, ylim=(-2,12))
plot!(years, collect(rf_high[1,:]), label="High Emissions", color=:purple, linewidth=2)
plot!(years, collect(rf_med[1,:]), label="Medium Emissions", color=:blue, linewidth=2)
plot!(years, collect(rf_low[1,:]), label="Low Emissions", color=:green, linewidth=2)

# combine all plots and format
all_plots3 = plot(f1, f2, layout=(1,2))
plot!(size=(1100,400), margin=5Plots.mm)
display(all_plots3)
savefig(all_plots3, "/Users/ced227/Desktop/plots/low_med_high_temp_rf.png")

# -------------------------------------------------------------------------------------------------------- #
# -------------- Plot 4: GMSLR Contributors for Low, Medium, High Emissions Comparison ------------------- #
# -------------------------------------------------------------------------------------------------------- #

# need for plotting
labels=["Antarctic" "GSIC" "Greenland" "LW Storage" "Thermal Expansion"]

# create a stacked area plot for sea level rise for LOW emissions
fig1 = plot(title="Low Emissions", xlabel="Year", ylabel="Global Mean Sea Level Anomaly (m)", legend=:topleft, xticks=xticks, ylim=(-1,10))
# create a matrix of sea level rise contributors
slr_contributions_low = [collect(antarctic_low[1,:])';
                         collect(gsic_low[1,:])';
                         collect(greenland_low[1,:])';
                         collect(lw_low[1,:])';
                         collect(te_low[1,:])']
areaplot!(years, slr_contributions_low', label=labels, color_palette=:darkrainbow, alpha=1, fillalpha=0.5)
plot!(years, collect(gmslr_low[1,:]), label="GMSLR", color=:black, linewidth=4, linestyle=:dash)

# create a stacked area plot for sea level rise for MEDIUM emissions
fig2 = plot(title="Medium Emissions", xlabel="Year", ylabel="Global Mean Sea Level Anomaly (m)", legend=:topleft, xticks=xticks, ylim=(-1,10))
# create a matrix of sea level rise contributors
slr_contributions_med = [collect(antarctic_med[1,:])';
                         collect(gsic_med[1,:])';
                         collect(greenland_med[1,:])';
                         collect(lw_med[1,:])';
                         collect(te_med[1,:])']
areaplot!(years, slr_contributions_med', label=labels, color_palette=:darkrainbow, alpha=1, fillalpha=0.5)
plot!(years, collect(gmslr_med[1,:]), label="GMSLR", color=:black, linewidth=4, linestyle=:dash)

# create a stacked area plot for sea level rise for HIGH emissions
fig3 = plot(title="High Emissions", xlabel="Year", ylabel="Global Mean Sea Level Anomaly (m)", legend=:topleft, xticks=xticks, ylim=(-1,10))
# create a matrix of sea level rise contributors
slr_contributions_high = [collect(antarctic_high[1,:])';
                         collect(gsic_high[1,:])';
                         collect(greenland_high[1,:])';
                         collect(lw_high[1,:])';
                         collect(te_high[1,:])']
areaplot!(years, slr_contributions_high', label=labels, color_palette=:darkrainbow, alpha=1, fillalpha=0.5)
plot!(years, collect(gmslr_high[1,:]), label="GMSLR", color=:black, linewidth=4, linestyle=:dash)

# combine all plots and format
all_plots4 = plot(fig1, fig2, fig3, layout=(1,3))
plot!(size=(1600,500), margin=10Plots.mm)
display(all_plots4)
savefig(all_plots4, "/Users/ced227/Desktop/plots/low_med_high_slr_contributions.png")

# -------------------------------------------------------------------------------------------------------- #
# --------------------- Supplementary: Can be useful for benchmarking ------------------------------------ #
# -------------------------------------------------------------------------------------------------------- #

"""# subset df to just include emissions parameters
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
upper_samples  = emission_params[upper_idx,:]"""