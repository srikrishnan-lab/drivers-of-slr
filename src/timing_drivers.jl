import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSVFiles
using Distributed
ncores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
print("Using $(ncores) cores...")
addprocs(div(ncores, 2))

@everywhere begin
    using DataFrames
    using Statistics
    using MLJ
    using EvoTrees
    using ShapML
end

EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees

# Assume that we call this script from the project folder
output_dir = "output/default"
slr_out = DataFrame(CSVFiles.load(joinpath(output_dir, "gmslr.csv")))
emissions = DataFrame(CSVFiles.load(joinpath(output_dir, "emissions.csv")))
parameters = DataFrame(CSVFiles.load(joinpath(output_dir, "parameters.csv")))

features = parameters

# make regression trees and compute Shapley values
function shapley_reg(levels, features)
    shap_df = DataFrame(feature_name = names(features))
    for lev in levels
        
        slr_reg_tree = EvoTreeRegressor(nrounds=200, max_depth=5)
        # set year of relevant sea level
        idx = findfirst.(>=(lev), eachrow(slr_out))
        reached_idx = (x -> !isnothing(x)).(idx) # find indices where sea level was reached
        idx = idx[reached_idx]
        params = features[reached_idx, :]
        yrs = parse.(Int, String.(idx))
        slr_reg_mach = machine(slr_reg_tree, params, yrs)
        MLJ.fit!(slr_reg_mach, force=true)
        # define function for parallelized Shapley calculation
        @everywhere function predict_slr(model, data)
            pred = DataFrame(slr_pred = MLJ.predict(model,data))
            return pred
        end

        explain = copy(features)
        reference = copy(features)
        shap_out = ShapML.shap(explain = explain, 
                                reference = reference,
                                model = slr_reg_mach,
                                predict_function = predict_slr,
                                sample_size = 200,
                                parallel = :samples,
                                seed = 1)
        shap_grouped = groupby(shap_out, :feature_name)
        shap_summary = combine(shap_grouped, 
            :shap_effect => (x -> mean(abs.(x))))
        rename!(shap_summary, Dict(:shap_effect_function => Symbol("mean_$(lev)")))
        shap_df = innerjoin(shap_df, shap_summary, on=:feature_name)
    end
    return shap_df
end

levels = 0.5:0.1:3.0
shap_df = shapley_reg(levels, features)
save("output/shapley/timing_shapley_indices.csv", shap_df)