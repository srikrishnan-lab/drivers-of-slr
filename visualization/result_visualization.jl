#= this script produces 2 plots to illustrate the outputs from different emissions trajectories:
Plot 1: shows the results for emissions, temp, rf, and GMSLR for ONE run and how it compares to the 95% credible interval of all runs (does not consider low,med,high emissions)
Plot 2: shows how low,med,high emissions translate to rf, temp, and GMSLR (shows median and 95% CI for low,med,high outputs) =#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots
using StatsPlots
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

num_samples = size(parameters, 1)
i = 1 #rand(1:num_samples) # random run to look at (between 1 and num_samples)

# create variables for years and xticks for plotting
years = parse.(Int64, names(emissions))
xticks = first(years):50:last(years)
# establish indices for historical data and future values
historical = findall((in)(first(years):2021), first(years):last(years)) # just historical data (1850-2021)
future = findall((in)(2022:2300), first(years):last(years)) # just future data (not historical)

# -------------------------------------------------------------------------------------------------------- #
# ---------------------- Plot 1: Input/Output Visualization for Individual Run --------------------------- #
# -------------------------------------------------------------------------------------------------------- #

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
#savefig(all_plots1, "/Users/ced227/Desktop/plots/run$i-$output_dir.png")

# -------------------------------------------------------------------------------------------------------- #
# ----- Plot 2: Representative Low, Medium, & High Emissions Trajectories with Output Visualization ------ #
# -------------------------------------------------------------------------------------------------------- #

# four panels of plot
p1 = plot(title="Selected CO₂ Emissions Trajectories", xlabel="Year", ylabel="Total CO₂ Emissions (GtCO₂/yr)") # p1 = emissions
p2 = plot(title="Radiative Forcing", xlabel="Year", ylabel="Global Radiative Forcing (W/m²)", ylim=(-2,12)) # p2 = radiative forcing
p3 = plot(title="Temperature", xlabel="Year", ylabel="Global Mean Temperature Anomaly (K)", ylim=(-1,9)) # p3 = temperature
p4 = plot(title="Global Mean Sea Level", xlabel="Year", ylabel="Global Mean Sea Level Anomaly (m)") # p4 = global mean sea level

# initialize values
run_name = ["high", "medium", "low"]
labels = ["High", "Medium", "Low"]
colors = [:purple, :blue, :green] # high, med, low emissions colors

for (i,run) in enumerate(run_name)
    # create variable for current run's label
    current_label = labels[i]

    # get results for run
    co2       = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "emissions.csv")))
    rf        = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "radiative_forcing.csv")))
    temp      = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "temperature.csv")))
    gmslr     = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "gmslr.csv")))
    antarctic = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "antarctic.csv")))
    gsic      = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "gsic.csv")))
    greenland = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "greenland.csv")))
    lw        = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "lw_storage.csv")))
    te        = DataFrame(load(joinpath(@__DIR__, "..", "results", "$run", "thermal_expansion.csv")))

    # Plot 1: add emissions curve for current run
    plot!(p1, years[future], collect(co2[1,:][future]), label="$current_label Emissions", color=colors[i], linewidth=3)

    # for Plots 2, 3, and 4: add median and credible interval
    uncertain_outputs = [(p2,rf), (p3,temp), (p4,gmslr)]
    for (plot,output) in uncertain_outputs 
        # create quantiles for current run and current output (95% credible interval & median)
        quantiles = mapslices(x -> quantile(x, [0.025, 0.5, 0.975]), Matrix(output), dims=1)
        # add credible interval
        plot!(plot, years[future], quantiles[1,:][future], fillrange=quantiles[3,:][future], fillalpha=0.5, alpha=0, color=colors[i], label="95% CI $current_label Emissions")
        # add median line
        plot!(plot, years[future], quantiles[2,:][future], color=colors[i], label="Median $current_label Emissions", linewidth=2)
    end

    # for all plots, add historical data
    all_outputs = [(p1,co2), (p2,rf), (p3,temp), (p4,gmslr)] # all plots
    if i == length(run_name) # if we're on the last run
        for (plot,output) in all_outputs
            # add historical data for all plots
            scatter!(plot, years[historical], collect(output[1,:][historical]), label="Historical Data", color=:black, markersize=2)
            if plot == p1 # emissions plot
                # add in extreme RCP scenario emissions (GtCO₂)
                rcp26, rcp85 = rcp_emissions()
                scatter!(p1, rcp85[:,1], rcp85[:,2], label="RCP 8.5", color=:black, markersize=3, shape=:rect)
                scatter!(p1, rcp26[:,1], rcp26[:,2], label="RCP 2.6", color=:black, markersize=4, shape=:utriangle)
            elseif plot == p2 # radiative forcing plot
                # add in vertical line for radiative forcing in 2100 (defines the RCP scenario)
                vline!([2100], color=:black, linestyle=:dash, label=:false)
            end
        end
    end
end

# combine plots and display
all_plots2 = plot(p1, p2, p3, p4, layout=4, size=(1000,800), margin=5Plots.mm, xticks=xticks, legend=:topleft, legendfontsize=7)
display(all_plots2)
#savefig(all_plots2, "/Users/ced227/Desktop/plots/low_med_high_output.pdf")