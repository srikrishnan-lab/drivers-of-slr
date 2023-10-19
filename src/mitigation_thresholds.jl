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

RegTree = @load EvoTreeRegressor pkg=EvoTrees

output_dir = "output/default"
slr_out = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))
# normalize relative to 2000
idx_2000 = findfirst(names(slr_out) .== "2000")
for row in axes(slr_out,1 )
    foreach(col -> slr_out[row, col] -= slr_out[row, idx_2000], axes(slr_out, 2))
end

milestones = [0.5, 1.0, 1.5, 2.0]
milestone = 2.0
idx = findfirst.(>=(milestone), eachrow(slr_out))
idx = replace(idx, nothing => Symbol(2300))
yrs = parse.(Int, String.(idx))

features = parameters
#refit tree using just t_peak, growth rate, and decline rate
key_params = [:t_peak, :gamma_g, :gamma_d]
slr_key_tree = RegTree(nrounds=200, max_depth=3)
slr_key_mach = machine(slr_key_tree, features[:, key_params], yrs)
MLJ.fit!(slr_key_mach)
# set up grid for predictions
growth_bds = minimum(features.gamma_g):0.001:maximum(features.gamma_g)
decline_bds = minimum(features.gamma_g):0.001:maximum(features.gamma_d)
coords = Iterators.product(growth_bds,decline_bds)

key_features = [[2050.0, grid[1], grid[2]] for grid in coords]
key_feature_df = DataFrame(mapreduce(permutedims, vcat, key_features), key_params)
predict_2050 = MLJ.predict(slr_key_mach, key_feature_df)

key_features = [[2060.0, grid[1], grid[2]] for grid in coords]
key_feature_df = DataFrame(mapreduce(permutedims, vcat, key_features), key_params)
predict_2060 = MLJ.predict(slr_key_mach, key_feature_df)

fig = Figure(resolution=(900, 500), fontsize=12, figure_padding=1)
axmain = Axis(fig[1,1], xlabel="Emissions Growth", ylabel="Emissions Decline")
axmain2 = Axis(fig[1,2], xlabel="Emissions Growth", ylabel="Emissions Decline")

fmap = Makie.contourf!(axmain, key_feature_df[:, 2], key_feature_df[:, 3], predict_2050, linewidth=2, levels=2050:25:2300)

Makie.contourf!(axmain2, key_feature_df[:, 2], key_feature_df[:, 3], predict_2060, linewidth=2, levels=2050:25:2300)

Colorbar(fig[2,1:2],fmap, vertical=false, flipaxis=true)


CairoMakie.save("figures/threshold_map.png", fig)
