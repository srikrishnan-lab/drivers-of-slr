# activate the environment
using Pkg
Pkg.activate(@__DIR__) # this does the same thing: Pkg.activate("/Users/ced227/Documents/Repositories/drivers-of-sea-level-rise/")
Pkg.instantiate()

# import necessary packages
using Mimi
using MimiDICE2013
using MimiSNEASY
using MimiBRICK
using DataFrames
using Impute
using CSVFiles
using Dates
using Distributions
using CSV

# include functions from other scripts
include("functions.jl")
include("create_dice_sneasy_brick.jl")

m = create_dice_sneasy_brick(start_year=2010, end_year=2305)
run(m)

# store dataframe with info for each RCP scenario
let 
    RCP26 = scenario_def("RCP26")
    RCP45 = scenario_def("RCP45")
    RCP60 = scenario_def("RCP60")
    RCP85 = scenario_def("RCP85")

#=
# write function so that normal random samples are positive
function pos_normal(μ,σ)
    x = rand(Normal(μ,σ))
    return (x ≥ 0) ? x : pos_normal(μ,σ)
end

run = @defsim begin 
    rv(co2_emissions) = pos_normal(μ=17.5, σ=8.4)
    rv(name2)
    rv(name3)
end
=#

for i in 1:5
    # run DICE
    dice = MimiDICE2013.get_model()
    run(dice)
    #explore(dice)

    # isolate output from DICE
    co2_emissions = getdataframe(dice, :emissions=>:E) # E is Total CO₂ emissions (GtCO₂ per year) (5 year timesteps)

    # create new df to fill in missing values for 1-year timesteps
    interpolated_df = DataFrame(Year = 2010:2305, Emissions = missing)
    # convert Emissions column type to both Missing and Float so that Float values can replace missing values
    interpolated_df[!, :Emissions] = convert(Vector{Union{Float64,Missing}}, interpolated_df[:, :Emissions])

    # loop through DICE df output
    for i in 1:nrow(co2_emissions)
        emission_value = co2_emissions[i, :E]
        loc = findfirst(==(co2_emissions[i, :time]), interpolated_df[:, :Year])
        interpolated_df[loc, :Emissions] = emission_value
    end

    # interpolate to fill in missing values in interpolated_df
    final_co2_df = Impute.interp(interpolated_df)
    final_emissions = final_co2_df[: , :Emissions]

    # create SNEASY-BRICK
    sneasy_brick = MimiBRICK.create_sneasy_brick(start_year=2010, end_year=2305)

    # update parameters to feed DICE emissions to SNEASY-BRICK
    update_param!(sneasy_brick, :ccm, :CO2_emissions, final_emissions)
    run(sneasy_brick)
    #explore(sneasy_brick)

    # retrieve output from SNEASY-BRICK
    slr_component = (:sea_level_rise, :slr_antartic_icesheet, :slr_glaciers_small_ice_caps, :slr_greeland_icesheet, :slr_landwater_storage, :slr_thermal_expansion)
    sea_level_rise = getdataframe(sneasy_brick, :global_sea_level=>slr_component) # total sea level rise from all components (m) (includes landwater storage for projection periods)

    rf_info = sneasy_brick[:radiativeforcing, :rf] # global radiative forcing (top of atmosphere) (W/m^2)
    temp_info = sneasy_brick[:ccm, :temp] # global mean temperature anomaly (K), relative to preindustrial

    # final output dataframe
    output_df = sea_level_rise
    # insert info on CO₂ emissions, radiative forcing, and global temperature
    insertcols!(output_df, 2, :co2_emissions => final_emissions, :radiative_forcing => rf_info, :temperature => temp_info)

    # match inputs to closest RCP scenario for each year
    n2o_concentration, aerosol_forcing, other_forcing = match_inputs(RCP26, RCP45, RCP60, RCP85, final_co2_df)

    # update parameters to reflect RCP matching
    update_param!(sneasy_brick, :rfco2, :N₂O, n2o_concentration)
    update_param!(sneasy_brick, :radiativeforcing, :rf_aerosol, aerosol_forcing)
    update_param!(sneasy_brick, :radiativeforcing, :rf_other, other_forcing)
    run(sneasy_brick)
    #explore(sneasy_brick)

    # export output to a .csv file
    CSV.write(joinpath(@__DIR__, "run$i.csv"), output_df, header=true)
    print(output_df)
end

end

#= create and view the netCDF file
create_netcdf(output_df=?)
view_nc = Dataset(joinpath(@__DIR__, "slr_output.nc"), "r")
close(view_nc)
=#