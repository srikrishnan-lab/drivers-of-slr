# activate the environment
using Pkg
Pkg.activate(joinpath(@__DIR__, "sea-level-rise-env/")) # this does the same thing: Pkg.activate("/Users/ced227/Documents/Repositories/drivers-of-sea-level-rise/sea-level-rise-env/")
Pkg.instantiate()

# import necessary packages
using Mimi
using MimiDICE2013
using MimiSNEASY
using MimiBRICK
using DataFrames
using Impute

# run DICE
DICE = MimiDICE2013.get_model()
run(DICE)
#explore(DICE)

# isolate output from DICE
co2_emissions = getdataframe(DICE, :emissions=>:E) # E is Total CO2 emissions (GtCO2 per year) (5 year timesteps)

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
SNEASY_BRICK = MimiBRICK.create_sneasy_brick(start_year=2010, end_year=2305)

# update parameters to feed DICE emissions to SNEASY-BRICK
update_param!(SNEASY_BRICK, :ccm, :CO2_emissions, final_emissions)
run(SNEASY_BRICK)
#explore(SNEASY_BRICK)

# retrieve output from SNEASY-BRICK
slr_component = (:sea_level_rise, :slr_antartic_icesheet, :slr_glaciers_small_ice_caps, :slr_greeland_icesheet, :slr_landwater_storage, :slr_thermal_expansion)
sea_level_rise = getdataframe(SNEASY_BRICK, :global_sea_level=>slr_component) # total sea level rise from all components (m) (includes landwater storage for projection periods)

rf_info = SNEASY_BRICK[:radiativeforcing, :rf] # global radiative forcing (top of atmosphere) (W/m^2)
temp_info = SNEASY_BRICK[:ccm, :temp] # global mean temperature anomaly (K), relative to preindustrial

# final output dataframe
output_df = sea_level_rise
# insert info on CO₂ emissions, radiative forcing, and global temperature
insertcols!(output_df, 2, :co2_emissions => final_emissions, :radiative_forcing => rf_info, :temperature => temp_info)

#CSV.write("/home/fs01/ced227/drivers-of-sea-level-rise/output.csv", output_df,  header = false)
#print(output_df)



# create a function that extracts and stores relevant RCP scenario information
function scenario_def(rcp_scenario)
    
    # load emissions and forcing data
    rcp_emissions      = DataFrame(load(joinpath(@__DIR__, "data", rcp_scenario*"_emissions.csv"), skiplines_begin=36))
    rcp_concentrations = DataFrame(load(joinpath(@__DIR__, "data", rcp_scenario*"_concentrations.csv"), skiplines_begin=37))
    rcp_forcing        = DataFrame(load(joinpath(@__DIR__, "data", rcp_scenario*"_midyear_radforcings.csv"), skiplines_begin=58))

    # index for years we want to consider
    year_idx = [246:1:541;] # rows 246-541
    years = rcp_emissions.YEARS[year_idx] # years 2010-2305

    # calculate CO₂ emissions
    rcp_co2_emissions = (rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)[year_idx] * 3.67 # multiply by 3.67 to convert GtC/yr to GtCO₂/yr

    # get RCP N₂O concentrations (used for CO₂ radiative forcing calculations)
    rcp_n2o_concentration = rcp_concentrations.N2O[year_idx] # units in ppb

    # calculate exogenous radiative forcings
    rcp_aerosol_forcing = (rcp_forcing.TOTAER_DIR_RF .+ rcp_forcing.CLOUD_TOT_RF)[year_idx] # units in W/m₂
    rcp_other_forcing   = (rcp_forcing.TOTAL_INCLVOLCANIC_RF .- rcp_forcing.CO2_RF .- rcp_forcing.TOTAER_DIR_RF .- rcp_forcing.CLOUD_TOT_RF)[year_idx] # units in W/m₂

    df = DataFrame(year = years, # years 2010-2305
                   rcp_co2_emissions = rcp_co2_emissions, # units in GtCO₂/yr
                   rcp_n2o_concentration = rcp_n2o_concentration, # units in ppb
                   rcp_aerosol_forcing = rcp_aerosol_forcing, # units in W/m₂
                   rcp_other_forcing = rcp_other_forcing) # units in W/m₂
    
    return df

end

# store dataframe for each RCP scenario
RCP26 = scenario_def("RCP26")
RCP45 = scenario_def("RCP45")
RCP60 = scenario_def("RCP60")
RCP85 = scenario_def("RCP85")

# create a function to match N₂O concentrations, aerosol forcing, and other forcing to closest RCP
function match_params()

    # initialize variables
    n2o_concentration = []
    aerosol_forcing = []
    other_forcing = []
    co2_emissions = [] # for reference (will not be using this parameter)

    for year in 1:nrow(final_co2_df)

        # isolate emission info for each scenario
        DICE_val = final_co2_df[year, :Emissions]
        RCP26_val = RCP26[year, :rcp_co2_emissions]
        RCP45_val = RCP45[year, :rcp_co2_emissions]
        RCP60_val = RCP60[year, :rcp_co2_emissions]
        RCP85_val = RCP85[year, :rcp_co2_emissions]

        # put four RCP emissions values in an array
        RCP_array = [RCP26_val, RCP45_val, RCP60_val, RCP85_val]
        loc = findmin(abs.(RCP_array .- DICE_val))[2] # find index of minimum distance

        # initialize variables
        co2_emission_val = 0
        n2o_concentration_val = 0
        aerosol_forcing_val = 0
        other_forcing_val = 0

        # assign values based on closest RCP Scenario
        if loc == 1
            # use RCP 2.6 values
            co2_emission_val = RCP26[year, :rcp_co2_emissions]
            n2o_concentration_val = RCP26[year, :rcp_n2o_concentration]
            aerosol_forcing_val = RCP26[year, :rcp_aerosol_forcing]
            other_forcing_val = RCP26[year, :rcp_other_forcing]
        elseif loc == 2
            # use RCP 4.5 values
            co2_emission_val = RCP45[year, :rcp_co2_emissions]
            n2o_concentration_val = RCP45[year, :rcp_n2o_concentration]
            aerosol_forcing_val = RCP45[year, :rcp_aerosol_forcing]
            other_forcing_val = RCP45[year, :rcp_other_forcing]
        elseif loc == 3
            # use RCP 6.0 values
            co2_emission_val = RCP60[year, :rcp_co2_emissions]
            n2o_concentration_val = RCP60[year, :rcp_n2o_concentration]
            aerosol_forcing_val = RCP60[year, :rcp_aerosol_forcing]
            other_forcing_val = RCP60[year, :rcp_other_forcing]
        elseif loc == 4
            # use RCP 8.5 values
            co2_emission_val = RCP85[year, :rcp_co2_emissions]
            n2o_concentration_val = RCP85[year, :rcp_n2o_concentration]
            aerosol_forcing_val = RCP85[year, :rcp_aerosol_forcing]
            other_forcing_val = RCP85[year, :rcp_other_forcing]
        end

        push!(CARBON, co2_emission_val)
        push!(n2o_concentration, n2o_concentration_val)
        push!(aerosol_forcing, aerosol_forcing_val)
        push!(other_forcing, other_forcing_val)

        #=
        println("n20 concen: ", n2o_concentration)
        println("aerosol: ", aerosol_forcing)
        println("other forcing: ", other_forcing)

        push!(n2o_concentration, n2o_concentration_val)
        push!(aerosol_forcing, aerosol_forcing_val)
        push!(other_forcing, other_forcing_val)

        println("DICE:    ", DICE_val)
        println("RCP 2.6: ", RCP26_val)
        println("RCP 4.5: ", RCP45_val)
        println("RCP 6.0: ", RCP60_val)
        println("RCP 8.5: ", RCP85_val)
        =#
    end

    return CARBON

end

# create plot with line for CO₂ emissions associated with each RCP scenario
plot(RCP26.year, RCP26.rcp_co2_emissions, xlabel="Year", ylabel="CO₂ emissions (GtCO₂/year)", label="RCP 2.6")
plot!(RCP45.year, RCP45.rcp_co2_emissions, label="RCP 4.5")
plot!(RCP60.year, RCP60.rcp_co2_emissions, label="RCP 6.0")
plot!(RCP85.year, RCP85.rcp_co2_emissions, label="RCP 8.5")

# add line for DICE scenario
years = final_co2_df[:,:Year]
emissions = final_co2_df[:,:Emissions]
plot!(years, emissions, label="DICE Scenario")

# add dots for closest RCP scenario to the DICE scenario in question
scatter!(years, match_params(), markersize=2, label="Closest RCP")
savefig(joinpath(@__DIR__, "param_plot.pdf"))

#=
update_param!(m, :ccm, :CO2_emissions, rcp_co2_emissions)
update_param!(m, :rfco2, :N₂O, rcp_n2o_concentration)
update_param!(m, :radiativeforcing, :rf_aerosol, rcp_aerosol_forcing)
update_param!(m, :radiativeforcing, :rf_other, rcp_other_forcing)
=#