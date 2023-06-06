# this script plots the results from the Global Sensitivity Analysis

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures
include("../src/functions.jl")

# initialization
output_dir = "default"
years = ["2100", "2200", "2300"]

# read in results
first_order = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "first_order.csv")))
first_CI    = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "first_CI.csv")))
total_order = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "total_order.csv")))
total_CI    = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "total_CI.csv")))

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

# ---------------------------- Plot stacked area plots for stratefied sensitivity analysis -------------------- #

# create df with sums of first order sensitivities for each group for each year
group_sums = combine(groupby(first_order, :group), names(first_order, Not([:group, :parameter])) .=> sum, renamecols=false)
group_sums[:,2:end] .= ifelse.(group_sums[:,2:end] .< 0, 0, group_sums[:,2:end]) # treat negative sensitivities as zero for each group

# sum of first order sensitivities for each year
totals_df = combine(group_sums, Not(:group) .=> sum; renamecols=false)

# define labels for plotting
group_names = unique(first_order.group)
group_sums = zeros(length(group_names), length(years))

# sum the first order sensitivites for each group
for (i, item) in enumerate(group_names)
    for (j, year) in enumerate(years)
        group_sums[i,j] = sum(first_order[first_order.group .== item, "$year"])
    end
end

sums_df = DataFrame(group_sums, collect(years))
sum(eachrow(sums_df))
insertcols!(sums_df, 1, :group => group_names)

# create a stacked area plot for contributors to first order sensitivity
fig1 = plot(title="Low Emissions", xlabel="Year", ylabel="First Order Sensitivity")
# create a matrix of sea level rise contributors
slr_contributions_low = [collect(antarctic_low[1,:])';
                         collect(gsic_low[1,:])';
                         collect(greenland_low[1,:])';
                         collect(lw_low[1,:])';
                         collect(te_low[1,:])']
areaplot!(years, slr_contributions_low', label=labels, color_palette=:darkrainbow, alpha=1, fillalpha=0.5)
plot!(years, collect(gmslr_low[1,:]), label="GMSLR", color=:black, linewidth=4, linestyle=:dash)

# ---------------------------- Plot first order sensitivities ------------------------------------------------- #
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
#savefig(p, "/Users/ced227/Desktop/plots/stratefied_gsa/late_first_order.pdf")

# ---------------------------- Plot total order sensitivities ------------------------------------------------- #
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
#savefig(p, "/Users/ced227/Desktop/plots/stratefied_gsa/late_total_order.pdf")