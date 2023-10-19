# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles # read CSV of Shapley indices
using DataFrames # data structure for indices
using Plots # plotting library
using Measures # adjust margins with explicit measures
using Statistics # get mean function

# assume this is called from the project root directory
output_path = "output/shapley"
shap_ind = DataFrame(CSVFiles.load(joinpath(output_path, "timing_shapley_indices.csv")))

# assign parameters to groups by subsystem
groups = Dict("Emissions"         => ["gamma_g", "t_peak", "gamma_d"],
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
group_col = Array{String}(undef, size(shap_ind)[1]) # preallocate space
for (key,value) in groups # loop through dictionary and create vector with group labels
    indices = findall((in)(value), shap_ind.feature_name) # find indices for current group
    group_col[indices] .= key # fill the indices with the name of current group
end
shap_ind.group = group_col # add group column to dataframe    

# average shapley effects by group for each year
shap_group = groupby(shap_ind[!, Not(:feature_name)], :group)
shap_group = combine(shap_group, Not(:group) .=> mean)
# normalize so the grouped shapley sums equal 1
shap_norm = mapcols(x -> x / sum(x), shap_group[!, Not([:group])])
insertcols!(shap_norm, 1, :group => shap_group.group)
shap_permute = permutedims(shap_norm, 1)

# plot indices over time
levels = 0.5:0.1:3.0
plt = areaplot(levels, Matrix(shap_permute[!, Not(:group)]), label=permutedims(shap_norm.group), xlabel="Sea Level Milestone (m)", ylabel="Normalized Grouped Shapley Index", color_palette=:seaborn_colorblind, left_margin=10mm, bottom_margin=10mm)
xticks!(plt, 0.5:0.5:3.0)
savefig(plt, "figures/timing-stacked-shapley-index.png")
