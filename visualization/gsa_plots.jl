# this script plots the results from the Global Sensitivity Analysis

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures
include("../src/functions.jl")

output_dir = "10_000_redo"

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
              "Statistical Noise"   => ["sd_temp", "sd_ocean_heat", "sd_glaciers", "sd_greenland", "sd_antarctic",
                                        "sd_gmsl", "rho_temperature", "rho_ocean_heat", "rho_glaciers", "rho_greenland",
                                        "rho_antarctic", "rho_gmsl"],
              "Climate System"      => ["CO2_0", "N2O_0", "temperature_0", "ocean_heat_0", "heat_diffusivity",
                                        "rf_scale_aerosol", "climate_sensitivity"],
              "Carbon Cycle"        => ["Q10", "CO2_fertilization", "CO2_diffusivity"],
              "Thermal Expansion"   => ["thermal_alpha", "thermal_s0"],
              "Greenland Ice Sheet" => ["greenland_a", "greenland_b", "greenland_alpha", "greenland_beta", "greenland_v0"],
              "Glaciers and SIC"    => ["glaciers_beta0", "glaciers_n", "glaciers_v0", "glaciers_s0"],
              "Antarctic Ocean"     => ["anto_alpha", "anto_beta"],
              "Antarctic Ice Sheet" => ["antarctic_s0", "antarctic_gamma", "antarctic_alpha", "antarctic_mu", "antarctic_nu",
                                        "antarctic_precip0", "antarctic_kappa", "antarctic_flow0", "antarctic_runoff_height0",
                                        "antarctic_c", "antarctic_bed_height0", "antarctic_slope", "antarctic_lambda",
                                        "antarctic_temp_threshold"],
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

# -------------------------------------- Plot first order sensitivities --------------------------------------- #
all_plts = []
years = ["2100", "2200", "2300"]

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
#savefig(p, "/Users/ced227/Desktop/plots/first_order.pdf")

# -------------------------------------- Plot total order sensitivities --------------------------------------- #
all_plts = []
years = ["2100", "2200", "2300"]

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
#savefig(p, "/Users/ced227/Desktop/plots/total_order.pdf")