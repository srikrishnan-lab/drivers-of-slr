#= this script divides the results of model_ensemble.jl into three buckets based on peaking time: before 2050, 2050-2100, and after 2100.
Then plots boxplots for relevant outputs for each of the three buckets for different years (either on centennial or decadal time scales). =#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures, StatsPlots
include("../src/functions.jl")

output_dir = "default"

# uncomment one of the two vectors below to select years
#years = ["2100", "2200", "2300"]
#years = ["2070", "2100", "2120", "2150"]

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

# establish indices for different buckets of peaking times
early_indices   = findall(parameters.t_peak .< 2050)           # before 2050
middle_indices  = findall(2050 .<= parameters.t_peak .< 2100)  # at least 2050 and before 2100
late_indices    = findall(parameters.t_peak .>= 2100)          # at least 2100

# --------------------------------------------------------------------------------- #
# -------------- Plot radiative forcing, temperature, and GMSLR ------------------- #
# --------------------------------------------------------------------------------- #
let
    all_outputs = []
    outputs = [radiative_forcing, temperature, gmslr]
    # loop through outputs of rf, temp, and GMSLR
    for (i, output) in enumerate(outputs)
        all_years = []
        # loop through years
        for year in years
            # boxplots for early, middle, and late peaking for current year
            titles = ["Radiative Forcing: $year" "Temperature: $year" "Global Mean Sea Level: $year"]
            ylabels = ["Global Radiative Forcing (W/m²)", "Global Mean Temperature Anomaly (K)", "Global Mean Sea Level Anomaly (m)"]
            xticks = (1:3, ["Early" "Middle" "Late"])
            ylims = []
            if length(years) == 3
                ylims = [(0, 15.5), (0, 11), (-0.5, 13.5)] # ylims for 2100, 2200, 2300
            elseif length(years) == 4
                ylims = [(0, 13), (0, 8.5), (-0.3, 5.5)] # ylims for 2070, 2100, 2120, 2150
            end
            p1 = boxplot(xticks=xticks, ylim=ylims[i], title=titles[i], ylabel=ylabels[i])
            boxplot!(output[early_indices,"$year"], color=:seagreen4) # early peaking
            boxplot!(output[middle_indices,"$year"], color=:lightblue) # middle peaking
            boxplot!(output[late_indices,"$year"], color=:plum4) # late peaking
            # add current year's plot to vector
            push!(all_years, p1)
        end 
        # combined plots for all years
        p2 = plot(all_years..., layout=(1,length(years)), legend=:false)
        push!(all_outputs, p2)
    end
    # set up values that change depending on number of years plotted
    title_size  = length(years)==3 ? 16 : 13 # title font size
    size        = length(years)==3 ? 12 : 11 # label/tick font size
    # combined plots for all outputs
    p3 = plot(all_outputs..., layout=(3,1), size=(1800,1200), bottom_margin=5mm, left_margin=10mm, ann=((0,1.06),:auto), dpi=600,
              titlefontsize=title_size, tickfontsize=size, labelfontsize=size)
    display(p3)
    #savefig(p3, "/Users/ced227/Desktop/plots/rf_temp_gmslr_boxplots.png")
end

# --------------------------------------------------------------------------------- #
# ------------------------- Plot contributors to GMSLR ---------------------------- #
# --------------------------------------------------------------------------------- #
let
    all_years = []
    buckets = [early_indices, middle_indices, late_indices]
    # loop through years
    for (i, year) in enumerate(years)
        all_buckets = []
        # loop through buckets of early, middle, and late peaking
        for (j, bucket) in enumerate(buckets)
            # boxplots for contributors to GMSLR
            titles = ["Early Peaking: $year" "Middle Peaking: $year" "Late Peaking: $year"]
            xticks = (1:5, ["Antarctic" "GSIC" "Greenland" "LW Storage" "TE"])
            ylims = []
            if length(years) == 3
                ylims = [(-0.4, 1.75), (-0.9, 5.3), (-1.2, 7.2)] # ylims for 2100, 2200, 2300
            elseif length(years) == 4
                ylims = [(-0.3, 0.9), (-0.4, 1.75), (-0.5, 2.4), (-0.8, 3.5)] # ylims for 2070, 2100, 2120, 2150
            end
            p1 = boxplot(xticks=xticks, ylims=ylims[i], title=titles[j], ylabel="Global Mean Sea Level Anomaly (m)")
            boxplot!(antarctic[bucket,"$year"])         # Antarctic ice sheet
            boxplot!(gsic[bucket,"$year"])              # glaciers and small ice caps
            boxplot!(greenland[bucket,"$year"])         # Greenland ice sheet
            boxplot!(lw_storage[bucket,"$year"])        # land water storage
            boxplot!(thermal_expansion[bucket,"$year"]) # thermal expansion
            # add current bucket's plot to vector
            push!(all_buckets, p1)
        end 
        # combined plots for early, middle, and late peaking
        p2 = plot(all_buckets..., layout=(1,3), legend=:false)
        push!(all_years, p2)
    end
    # set up values that change depending on number of years plotted
    loc         = length(years)==3 ? (0,1.06) : (0,1.09) # location of panel labels
    title_size  = length(years)==3 ? 16 : 15 # title font size
    size        = length(years)==3 ? 11 : 10 # label/tick font size
    # combined plots for all years
    p3 = plot(all_years..., layout=(length(years),1), size=(1800,1200), bottom_margin=5mm, left_margin=10mm, ann=(loc,:auto), dpi=600,
              titlefontsize=title_size, tickfontsize=size, labelfontsize=size)
    display(p3)
    #savefig(p3, "/Users/ced227/Desktop/plots/slr_boxplots.png")
end

# -------------------------------------------------------------------------------------------------------- #
# -------------- Analysis of Parameter/Outcome Percentiles for each Peaking Group ------------------------ #
# -------------------------------------------------------------------------------------------------------- #
let
    # intialize variables and df
    groups = [early_indices, middle_indices, late_indices]
    group_names = ["early", "middle", "late"]
    param_names = names(parameters)
    all_percentiles = DataFrame([group => zeros(length(param_names)) for group in group_names])
    insertcols!(all_percentiles, 1, :parameter => param_names)

    # different percentiles to consider
    percentile_emissions    = zeros(length(groups))
    percentile_temps        = zeros(length(groups))
    percentile_gmslr        = zeros(length(groups))
    percentile_greenland    = zeros(length(groups))
    percentile_antarctic    = zeros(length(groups))
    percentile_gsic         = zeros(length(groups))
    percentile_lws          = zeros(length(groups))
    percentile_te           = zeros(length(groups))

    for (i, group) in enumerate(groups) # loop through peaking group indices
        outcome = gmslr # change desired outcome as needed (i.e., gmslr, antarctic)
        # establish index for maximum value of outcome in 2300 for current group
        idx = findfirst(maximum(outcome[group, end]) .∈ outcome[group, end])
        # print maximum value for each group
        current_group = group_names[i]
        println("$current_group: ", maximum(outcome[group, end]))
        # print cumulative emissions for each group
        cum_emissions = sum.(eachrow(emissions[group,:])) # cumulative emissions in 2300 for current peaking group
        println("$current_group-->emissions: ", cum_emissions[idx])
        
        for (j, param) in enumerate(param_names) # loop through parameters
            # find the value and percentile for the sample
            params_subset = parameters[group,:] # df just with current peaking group
            value = params_subset[idx, param]
            percentile = 100 * mean(params_subset[:, param] .< value)
            # add the percentile to df
            all_percentiles[j,i+1] = percentile
        end

        # add percentile for each metric in 2300
        percentile_emissions[i] = 100 * mean(cum_emissions .< cum_emissions[idx])
        percentile_temps[i]     = 100 * mean(temperature[group,end] .< temperature[group,end][idx])
        percentile_gmslr[i]     = 100 * mean(gmslr[group,end] .< gmslr[group,end][idx])
        percentile_greenland[i] = 100 * mean(greenland[group,end] .< greenland[group,end][idx])
        percentile_antarctic[i] = 100 * mean(antarctic[group,end] .< antarctic[group,end][idx])
        percentile_gsic[i]      = 100 * mean(gsic[group,end] .< gsic[group,end][idx])
        percentile_lws[i]       = 100 * mean(lw_storage[group,end] .< lw_storage[group,end][idx])
        percentile_te[i]        = 100 * mean(thermal_expansion[group,end] .< thermal_expansion[group,end][idx])

        if i == length(groups) # if we're on the last run
            # insert each percentile vector into first row of all_percentiles df 
            insert!.(eachcol(all_percentiles), 1, ["cumulative_emissions", percentile_emissions...])
            insert!.(eachcol(all_percentiles), 1, ["temperature", percentile_temps...])
            insert!.(eachcol(all_percentiles), 1, ["gmslr", percentile_gmslr...])
            insert!.(eachcol(all_percentiles), 1, ["greenland_slr", percentile_greenland...])
            insert!.(eachcol(all_percentiles), 1, ["antarctic_slr", percentile_antarctic...])
            insert!.(eachcol(all_percentiles), 1, ["gsic_slr", percentile_gsic...])
            insert!.(eachcol(all_percentiles), 1, ["lw_storage_slr", percentile_lws...])
            insert!.(eachcol(all_percentiles), 1, ["thermal_expansion_slr", percentile_te...])
        end
    end
    show(all_percentiles, allrows=true)
end