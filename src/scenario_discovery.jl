# scenario_discovery.jl: code to take SLR ensemble and construct a classification tree to find high/low outcomes

import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles
using Statistics
using Plots
using DataFrames
using MLJ
using EvoTrees
using ShapML
using CairoMakie


EvoTreeClassifier = @load EvoTreeClassifier pkg=EvoTrees

# Assume that we call this script from the project folder
output_dir = "output/default"
slr_out = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
# normalize relative to 2000
idx_2000 = findfirst(names(slr_out) .== "2000")
for row in axes(slr_out,1 )
    foreach(col -> slr_out[row, col] -= slr_out[row, idx_2000], axes(slr_out, 2))
end


## Start with SLR outcomes in 2100
yr = 2100

features = parameters

# Label SLR outcomes as low and high based on quantiles
hi_threshold = 1.0
slr_hi_labels = ifelse.(hi_threshold .< slr_out[:, Symbol(yr)], "high", "normal")

# classify extreme high outcomes
slr_hi_tree = EvoTreeClassifier(nrounds=200, max_depth=3)
slr_hi_mach = machine(slr_hi_tree, features, slr_hi_labels)
MLJ.fit!(slr_hi_mach, force=true)

imp_out = stack(DataFrame(feature_importances(slr_hi_mach)))
fig_imp2100 = Figure()
ax = Axis(fig_imp2100[1,1], xticks = (1:10, imp_out.variable[1:10]), xticklabelrotation = pi/4, ylabel = "Normalized Shapley Index")
imp_plt = Makie.barplot!(ax, imp_out.value[1:10])
save("figures/importance_bar_2100.png", fig_imp2100)

#refit tree using just t_peak, antarctic temp threshold, and climate sensitivity
key_params = [:t_peak, :antarctic_temp_threshold, :climate_sensitivity]
slr_key_tree = EvoTreeClassifier(nrounds=200, max_depth=3)
slr_key_mach = machine(slr_key_tree, features[:, key_params], slr_hi_labels)
MLJ.fit!(slr_key_mach)
# set up grid for predictions
temp_bds = round(minimum(features.antarctic_temp_threshold), digits=1):0.001:round(maximum(features.antarctic_temp_threshold), digits=1)
ecs_bds = round(minimum(features.climate_sensitivity), digits=1):0.001:round(maximum(features.climate_sensitivity),digits=1)
coords = Iterators.product(temp_bds,ecs_bds)

# plot dividing contour
key_features = [[2050.0, grid[1], grid[2]] for grid in coords]
key_feature_df = DataFrame(mapreduce(permutedims, vcat, key_features), key_params)
predict_key = MLJ.predict(slr_key_mach, key_feature_df)
predict_class_2050 = [pred.prob_given_ref[1] for pred in predict_key]

key_features = [[2070.0, grid[1], grid[2]] for grid in coords]
key_feature_df = DataFrame(mapreduce(permutedims, vcat, key_features), key_params)
predict_key = MLJ.predict(slr_key_mach, key_feature_df)
predict_class_2070 = Float64.([pred.prob_given_ref[1] for pred in predict_key])

key_features = [[2090.0, grid[1], grid[2]] for grid in coords]
key_feature_df = DataFrame(mapreduce(permutedims, vcat, key_features), key_params)
predict_key = MLJ.predict(slr_key_mach, key_feature_df)
predict_class_2090 = Float64.([pred.prob_given_ref[1] for pred in predict_key])

contour_colors = cgrad(:Reds_5)

fig = Figure(resolution= (350, 300), fontsize=12, figure_padding=1)
axmain = Axis(fig[2,1], xlabel="Equilibrium Climate Sensitivity (°C)", ylabel="Antarctic Temperature Threshold (°C)")
axtop = Axis(fig[1,1], limits=((1.5, 6), (0, nothing)))
axright = Axis(fig[2, 2], limits=((0, nothing), (1.3, 3.5)))
linkyaxes!(axmain, axright)
linkxaxes!(axmain, axtop)

fmap = Makie.contour!(axmain, key_feature_df[:, 3], 15.42 .+ 0.8365 *key_feature_df[:, 2], predict_class_2050, color=contour_colors[0.2], levels=0.0:0.5:1.0, linewidth=4, xlabel="Equilibrium Climate Sensitivity (°C)", ylabel="Antarctic Temperature Threshold (°C)")
fmap2 =Makie.contour!(axmain, key_feature_df[:, 3], 15.42 .+ 0.8365 *key_feature_df[:, 2], predict_class_2070, color=contour_colors[0.5], levels=0.0:0.5:1.0, linewidth=4)
fmap3 = Makie.contour!(axmain, key_feature_df[:, 3], 15.42 .+ 0.8365 *key_feature_df[:, 2], predict_class_2090, color=contour_colors[0.8], levels=0.0:0.5:1.5, linewidth=4)

Makie.ylims!(1.25, 3.75)

text!(axmain, 3.25, 1.4, text=">50% Probability of \n GMSLR >1m in 2100")
text!(axmain, 2.0, 3.1, text="<50% Probability of \n GMSLR >1m in 2100")

# plot marginal densities for the geophysical uncertainties
Makie.density!(axtop, parameters.climate_sensitivity)
Makie.density!(axright, 15.42 .+ 0.8365 * parameters.antarctic_temp_threshold,direction=:y)
hidedecorations!(axtop)
hidedecorations!(axright)
hidespines!(axtop)
hidespines!(axright)

# create legend
elem_2050 = LineElement(color = contour_colors[0.2], linestyle = nothing, linewidth=4)
elem_2070 = LineElement(color = contour_colors[0.5], linestyle = nothing, linewidth=4)
elem_2090 = LineElement(color = contour_colors[0.8], linestyle = nothing, linewidth=4)
leg = Legend(fig[3,1], [elem_2050, elem_2070, elem_2090], ["2050", "2070", "2090"], "Year Emissions Peak", orientation=:horizontal, tellwidth=false, tellheight=true, framevisible=false)

colsize!(fig.layout, 2, Auto(0.25))
rowsize!(fig.layout, 1, Auto(0.25))
colgap!(fig.layout, 1, Relative(0.0))
rowgap!(fig.layout, 1, Relative(0.0))
rowgap!(fig.layout, 2, Relative(0.01))

CairoMakie.save("figures/factor_map.png", fig)


