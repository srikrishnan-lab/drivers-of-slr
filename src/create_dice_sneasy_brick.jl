# include functions from other script
include("functions.jl")

# function returns a Mimi model instance with MimiDICE, MimiSNEASY, and MimiBRICK linked together
function create_dice_sneasy_brick(; start_year::Integer=2010, end_year::Integer=2305, DICE_model=DICE, SNEASY_BRICK_model=SNEASY_BRICK)

    # ------------------------------- Generate emissions from DICE and feed them to SNEASY-BRICK ------------------------------------ #

    # isolate output from DICE
    co2_emissions = getdataframe(DICE, :emissions=>:E) # E is Total CO₂ emissions (GtCO₂ per year) (5 year timesteps)

    # create new df to fill in missing values for 1-year timesteps
    interpolated_df = DataFrame(Year = 2010:2305, Emissions = missing)
    # convert Emissions column type to both Missing and Float so that Float values can replace missing values
    interpolated_df[!, :Emissions] = convert(Vector{Union{Float64,Missing}}, interpolated_df[:, :Emissions])

    # loop through DICE df output, results in df with a value every 5 years, Missing in other years
    for i in 1:nrow(co2_emissions)
        emission_value = co2_emissions[i, :E]
        loc = findfirst(==(co2_emissions[i, :time]), interpolated_df[:, :Year])
        interpolated_df[loc, :Emissions] = emission_value
    end

    # interpolate to fill in Missing values in interpolated_df
    final_co2_df = Impute.interp(interpolated_df)
    final_emissions = final_co2_df[: , :Emissions]

    # update parameters to feed DICE emissions to SNEASY-BRICK
    update_param!(SNEASY_BRICK, :ccm, :CO2_emissions, final_emissions)

    # -------------- Match inputs for N₂O concentration, aerosol forcing, and other forcing to closest RCP scenario ----------------- #

    # store dataframe with info for each RCP scenario
    RCP26 = scenario_def("RCP26")
    RCP45 = scenario_def("RCP45")
    RCP60 = scenario_def("RCP60")
    RCP85 = scenario_def("RCP85")

    # match inputs to closest RCP scenario for each year
    n2o_concentration, aerosol_forcing, other_forcing = match_inputs(RCP26, RCP45, RCP60, RCP85, final_co2_df)

    # update inputs to reflect RCP matching
    update_param!(SNEASY_BRICK, :rfco2, :N₂O, n2o_concentration)
    update_param!(SNEASY_BRICK, :radiativeforcing, :rf_aerosol, aerosol_forcing)
    update_param!(SNEASY_BRICK, :radiativeforcing, :rf_other, other_forcing)

    # return model ensemble
    return DICE, SNEASY_BRICK
end




#----------------------------I am using the below for testing, can delete later----------------------------#
using MimiDICE2013
DICE = MimiDICE2013.get_model()
run(DICE)

# isolate output from DICE
co2_emissions = getdataframe(DICE, :emissions=>:E) # E is Total CO₂ emissions (GtCO₂ per year) (5 year timesteps)

# create new df to fill in missing values for 1-year timesteps
interpolated_df = DataFrame(Year = 2010:2305, Emissions = missing)
# convert Emissions column type to both Missing and Float so that Float values can replace missing values
interpolated_df[!, :Emissions] = convert(Vector{Union{Float64,Missing}}, interpolated_df[:, :Emissions])

# loop through DICE df output, results in df with a value every 5 years, Missing in other years
for i in 1:nrow(co2_emissions)
    emission_value = co2_emissions[i, :E]
    loc = findfirst(==(co2_emissions[i, :time]), interpolated_df[:, :Year])
    interpolated_df[loc, :Emissions] = emission_value
end

# interpolate to fill in Missing values in interpolated_df
final_co2_df = Impute.interp(interpolated_df)
final_emissions = final_co2_df[: , :Emissions]

# update parameters to feed DICE emissions to SNEASY-BRICK
update_param!(SNEASY_BRICK, :ccm, :CO2_emissions, final_emissions)

# final_emissions goes from 2010-2305
length(final_emissions)