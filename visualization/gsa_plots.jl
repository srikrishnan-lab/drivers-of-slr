#= this script creates 3 plots and 1 .xlsx file showing the results from the Global Sensitivity Analysis:
Plot 1: 4-panel stacked area plot (one for each run) showing how the contributors to first order sensitivity change over time
Plot 2: 3-panel scatterplot (showing 3 years) that shows the first order sensitivity indices for each parameter for a specified run 
Plot 3: same as Plot 2, except for total order sensitivities
The last result is an .xlsx file with tables that show selected parameter interactions for 3 years: 2100, 2200, and 2300 =#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures
include("../src/functions.jl")

# initialization
run_name = "late" # choose from "default", "early", "middle", or "late"
years = ["2100", "2200", "2300"]

# read in first and total order results
first_order  = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$run_name", "first_order.csv")))
first_CI     = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$run_name", "first_CI.csv")))
total_order  = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$run_name", "total_order.csv")))
total_CI     = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$run_name", "total_CI.csv")))

# treat negative sensitivities as zero
first_order[:,2:end] .= ifelse.(first_order[:,2:end] .< 0, 0, first_order[:,2:end])
total_order[:,2:end] .= ifelse.(total_order[:,2:end] .< 0, 0, total_order[:,2:end])

# -------------------------- Group parameters and create column for plot color scheme --------------------------- #

# create groups for color scheme in plots
groups = Dict("Emissions"           => ["gamma_g", "t_peak", "gamma_d"],
              "Statistical Noise"   => ["sd_temp", "sd_ocean_heat", "sd_glaciers", "sd_greenland", "sd_antarctic", "sd_gmsl",
                                        "rho_temperature", "rho_ocean_heat", "rho_glaciers", "rho_greenland", "rho_antarctic", "rho_gmsl"],
              "Climate System"      => ["temperature_0", "ocean_heat_0", "heat_diffusivity", "rf_scale_aerosol", "climate_sensitivity"],
              "Carbon Cycle"        => ["CO2_0", "N2O_0", "Q10", "CO2_fertilization", "CO2_diffusivity"],
              "Thermal Expansion"   => ["thermal_alpha", "thermal_s0"],
              "Greenland"           => ["greenland_a", "greenland_b", "greenland_alpha", "greenland_beta", "greenland_v0"],
              "Glaciers and SIC"    => ["glaciers_beta0", "glaciers_n", "glaciers_v0", "glaciers_s0"],
              "Antarctic"           => ["antarctic_s0", "antarctic_gamma", "antarctic_alpha", "antarctic_mu", "antarctic_nu",
                                        "antarctic_precip0", "antarctic_kappa", "antarctic_flow0", "antarctic_runoff_height0",
                                        "antarctic_c", "antarctic_bed_height0", "antarctic_slope", "antarctic_lambda",
                                        "antarctic_temp_threshold", "anto_alpha", "anto_beta"],
              "Land Water Storage"  => ["lw_random_sample"])

# preallocate space for new column of parameter groupings
groups_col = Array{String}(undef, length(first_order.parameter))

# loop through each key-value pair in dictionary
for (key,value) in groups
    # find indices for current parameter group
    indices = findall((in)(value), first_order.parameter)
    # fill the indices with the name of the group
    groups_col[indices] .= key
end

# add a column for the parameter groupings
insertcols!(first_order, 1, :group => groups_col)
insertcols!(total_order, 1, :group => groups_col)

# ---------------------------- Plot first order sensitivities ------------------------------------------------- #
# plots for the run specified in "run_name"
all_plts = []

for year in years
    plt = scatter(first_order.parameter, first_order[:,"$year"], yerror = first_CI[:,"$year"], 
                  group=first_order.group, palette=:tab10, xticks=:all, xrotation=45,
                  title = "$year", xlabel="Parameter", ylabel="Sensitivity", legend=:true,
                  markersize=5, markerstrokewidth=0.5, tickfontsize=6, legendfontsize=6,
                  bottom_margin=7mm, left_margin=5mm)
    push!(all_plts, plt)
end
#ylim=(-0.03,0.82), yticks=0:0.1:0.8,

p = plot(all_plts..., plot_title="First Order Sensitivities", layout=(3,1), size=(900,1200))
display(p)
#savefig(p, "/Users/ced227/Desktop/plots/stratefied_gsa/UPDATE_first_order.pdf")

# ---------------------------- Plot total order sensitivities ------------------------------------------------- #
# plots for the run specified in "run_name"
all_plts = []

for year in years
    plt = scatter(total_order.parameter, total_order[:,"$year"], yerror = total_CI[:,"$year"], 
                  group=total_order.group, palette=:tab10, xticks=:all, xrotation=45,
                  title = "$year", xlabel="Parameter", ylabel="Sensitivity", legend=:true,
                  markersize=5, markerstrokewidth=0.5, tickfontsize=6, legendfontsize=6,
                  bottom_margin=7mm, left_margin=5mm)
    push!(all_plts, plt)
end
#ylim=(-0.02,0.62), yticks=0:0.1:0.6,

p = plot(all_plts..., plot_title="Total Order Sensitivities", layout=(3,1), size=(900,1200))
display(p)
#savefig(p, "/Users/ced227/Desktop/plots/stratefied_gsa/UPDATE_total_order.pdf")

# ---------------------------- Plot stacked area plots for stratefied sensitivity analysis -------------------- #
# plots for all runs: default, early, middle, and late
all_plts = []
gsa_runs = ["default", "early", "middle", "late"]
titles = ["Full Ensemble", "Early Peaking", "Middle Peaking", "Late Peaking"]

for (i, item) in enumerate(gsa_runs) # loop through each run
    # read in first order results for current run
    current_run = gsa_runs[i]
    gsa_first = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$current_run", "first_order.csv")))
    
    # add a column for parameter groups
    group_col = Array{String}(undef, length(gsa_first.parameter)) # preallocate space
    for (key,value) in groups # loop through dictionary
        indices = findall((in)(value), gsa_first.parameter) # find indices for current group
        group_col[indices] .= key # fill the indices with the name of current group
    end
    insertcols!(gsa_first, 1, :group => group_col) # add the column

    # create df with sums of first order sensitivities for each group for each year
    group_sums = combine(groupby(gsa_first, :group), names(gsa_first, Not([:group, :parameter])) .=> sum, renamecols=false)
    group_sums[:,2:end] .= ifelse.(group_sums[:,2:end] .< 0, 0, group_sums[:,2:end]) # treat negative sensitivities as zero for each group
    group_sums = group_sums[[8,3,4,1,7,6,9,2,5], :] # reorder rows for consistent color scheme

    # normalize sensitivities so that they add up to 1 each year
    totals_df = combine(group_sums, Not(:group) .=> sum; renamecols=false) # sum of first order sensitivities for each year
    normalized_df = group_sums[:,2:end] ./ totals_df # divide each group's value by the sum for each year
    insertcols!(normalized_df, 1, :group => group_sums.group) # re-add group column

    # intialize years and legend labels
    all_years = parse.(Int64, names(normalized_df)[2:end]) # years (x-axis)
    labels = permutedims(vcat(normalized_df.group)) # group names

    # create a stacked area plot for contributors to first order sensitivity
    plt = plot(title=titles[i], xlabel="Year", ylabel="First Order Index", bottom_margin=7mm, left_margin=5mm, legend=:right)
    areaplot!(all_years, Matrix(normalized_df[:,2:end])', label=labels, palette=:tab10, alpha=1, fillalpha=0.7, legendfontsize=6)
    push!(all_plts, plt) # add current plot to array
end

p = plot(all_plts..., plot_title="First Order Sensitivity Stacked Area Plots", layout=(2,2), size=(1200,900))
display(p)
#savefig(p, "/Users/ced227/Desktop/plots/stacked_area.pdf")

# ---------------------------- Tables of second order sensitivities -------------------------------------------- #
# plots for the "default" run

param_pairs = [("gamma_g","t_peak"), # vector of tuples with parameter pairs we want to consider
            ("t_peak","gamma_d"),
            ("gamma_d","sd_glaciers"), 
            ("gamma_d","rho_ocean_heat"),
            ("gamma_d","rho_gmsl"),
            ("gamma_d","greenland_v0"),
            ("gamma_d","Q10"),
            ("gamma_d","CO2_fertilization"),
            ("gamma_d","CO2_diffusivity"),
            ("gamma_d","anto_alpha"),
            ("gamma_d","anto_beta"),
            ("gamma_d","antarctic_c"),
            ("gamma_d","antarctic_mu"),
            ("gamma_d","antarctic_precip0"),
            ("gamma_d","antarctic_runoff_height0"),
            ("gamma_d","antarctic_bed_height0"),
            ("gamma_d","antarctic_lambda"),
            ("gamma_d","antarctic_alpha"),
            ("gamma_d","lw_random_sample")]

# preallocate storage for sensitivity indices and confidence interval bounds
estimates    = Matrix{Any}(undef, length(param_pairs), length(years))
lower_bounds = Matrix{Any}(undef, length(param_pairs), length(years))
upper_bounds = Matrix{Any}(undef, length(param_pairs), length(years))

for (i, year) in enumerate(years) # loop through years considered
    # get results for current year
    second_order = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "default", "second_order", "sens_$year.csv")))
    second_CI    = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "default", "second_order", "CI_$year.csv")))
    # treat negative sensitivities as zero
    second_order[:,2:end] .= ifelse.(second_order[:,2:end] .< 0, 0, second_order[:,2:end])

    for (j, pair) in enumerate(param_pairs) # loop through parameter pairs
        # separate pair into two parameters
        param1 = pair[1]
        param2 = pair[2]
        # get second order index
        gsa_value = second_order[(second_order.parameter .== "$param1"), "$param2"][1]
        # get the 95% confidence interval
        margin_of_error = second_CI[(second_CI.parameter .== "$param1"), "$param2"][1] # get margin of error from results
        lower_bound = gsa_value - margin_of_error
        upper_bound = gsa_value + margin_of_error
        # add in second order index and lower/upper bounds to matrices
        estimates[j,i]    = gsa_value
        lower_bounds[j,i] = lower_bound
        upper_bounds[j,i] = upper_bound
    end
end

# convert each Matrix to DataFrame
estimates_df    = DataFrame(estimates, years)
lower_bounds_df = DataFrame(lower_bounds, years)
upper_bounds_df = DataFrame(upper_bounds, years)

# add columns to each df for clarity
dfs = [estimates_df, lower_bounds_df, upper_bounds_df]
for df in dfs
    insertcols!(df, 1, :"Parameter 1" => first.(param_pairs)) # parameter 1 column
    insertcols!(df, 2, :"Parameter 2" => last.(param_pairs)) # parameter 2 column
end

# save to Excel workbook, with each df as a separate sheet
#XLSX.writetable("/Users/ced227/Desktop/plots/second_order.xlsx", "estimates" => estimates_df, "lower_bounds" => lower_bounds_df, "upper_bounds" => upper_bounds_df)