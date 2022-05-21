using CSVFiles
using DataFrames
using Mimi
using Plots
using MimiSNEASY

#--------------------------------------------------Plot 1--------------------------------------------------#

function plot_scenarios(rcp_scenario)

    # Load emissions and forcing data (index into appropriate years).
    rcp_emissions      = DataFrame(load(joinpath(@__DIR__, "data", rcp_scenario*"_emissions.csv"), skiplines_begin=36))
    rcp_concentrations = DataFrame(load(joinpath(@__DIR__, "data", rcp_scenario*"_concentrations.csv"), skiplines_begin=37))
    rcp_forcing        = DataFrame(load(joinpath(@__DIR__, "data", rcp_scenario*"_midyear_radforcings.csv"), skiplines_begin=58))

    # Calculate CO₂ emissions.
    rcp_co2_emissions = (rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)[[1:1:736;]]

    # Get RCP N₂O concentrations (used for CO₂ radiative forcing calculations).
    rcp_n2o_concentration = rcp_concentrations.N2O[[1:1:736;]]

    # Calculate exogenous radiative forcings.
    rcp_aerosol_forcing    = (rcp_forcing.TOTAER_DIR_RF .+ rcp_forcing.CLOUD_TOT_RF)[[1:1:736;]]
    rcp_other_forcing      = (rcp_forcing.TOTAL_INCLVOLCANIC_RF .- rcp_forcing.CO2_RF .- rcp_forcing.TOTAER_DIR_RF .- rcp_forcing.CLOUD_TOT_RF)[[1:1:736;]]

    df = DataFrame(year = [1765:1:2500;],
                   rcp_co2_emissions = rcp_co2_emissions, 
                   rcp_n2o_concentration = rcp_n2o_concentration,
                   rcp_aerosol_forcing = rcp_aerosol_forcing,
                   rcp_other_forcing = rcp_other_forcing)

    return df

end

# create dataframe for each RCP scenario
df_RCP26 = plot_scenarios("RCP26")
df_RCP45 = plot_scenarios("RCP45")
df_RCP60 = plot_scenarios("RCP60")
df_RCP85 = plot_scenarios("RCP85")

# years and labels for plots
years = [1765:1:2500;]
labels = ["RCP 2.6" "RCP 4.5" "RCP 6.0" "RCP 8.5"]

# plot rcp_co2_emissions
all_rcp_co2_emissions = [df_RCP26[:,:rcp_co2_emissions] df_RCP45[:,:rcp_co2_emissions] df_RCP60[:,:rcp_co2_emissions] df_RCP85[:,:rcp_co2_emissions]]
emissions_plot = plot(years, all_rcp_co2_emissions, label = labels)
plot!(emissions_plot, title = "rcp_co2_emissions", xlabel = "year", ylabel = "rcp_co2_emissions")

# plot rcp_n2o_concentration
all_rcp_n2o_concentration = [df_RCP26[:,:rcp_n2o_concentration] df_RCP45[:,:rcp_n2o_concentration] df_RCP60[:,:rcp_n2o_concentration] df_RCP85[:,:rcp_n2o_concentration]]
n2o_plot = plot(years, all_rcp_n2o_concentration, label = labels, legend = :topleft)
plot!(n2o_plot, title = "rcp_n2o_concentration", xlabel = "year", ylabel = "rcp_n2o_concentration")

# plot rcp_aerosol_forcing
all_rcp_aerosol_forcing = [df_RCP26[:,:rcp_aerosol_forcing] df_RCP45[:,:rcp_aerosol_forcing] df_RCP60[:,:rcp_aerosol_forcing] df_RCP85[:,:rcp_aerosol_forcing]]
aerosol_plot = plot(years, all_rcp_aerosol_forcing, label = labels)
plot!(aerosol_plot, title = "rcp_aerosol_forcing", xlabel = "year", ylabel = "rcp_aerosol_forcing")

# plot rcp_other_forcing
all_rcp_other_forcing = [df_RCP26[:,:rcp_other_forcing] df_RCP45[:,:rcp_other_forcing] df_RCP60[:,:rcp_other_forcing] df_RCP85[:,:rcp_other_forcing]]
other_forcing_plot = plot(years, all_rcp_other_forcing, label = labels, legend = :topleft)
plot!(other_forcing_plot, title = "rcp_other_forcing", xlabel = "year", ylabel = "rcp_other_forcing")

# display all plots
#=
display(emissions_plot)
display(n2o_plot)
display(aerosol_plot)
display(other_forcing_plot)
=#
all_plots = plot(emissions_plot, n2o_plot, aerosol_plot, other_forcing_plot, size=(1300,1100))
display(all_plots)

#--------------------------------------------------Plot 2--------------------------------------------------#

# Scenario 1: update all parameters to RCP 8.5
SNEASY = MimiSNEASY.get_model()
update_param!(SNEASY, :ccm, :CO2_emissions, df_RCP85[:, :rcp_co2_emissions])
update_param!(SNEASY, :rfco2, :N₂O, df_RCP85[:, :rcp_n2o_concentration]) 
update_param!(SNEASY, :radiativeforcing, :rf_aerosol, df_RCP85[:, :rcp_aerosol_forcing])
update_param!(SNEASY, :radiativeforcing, :rf_other, df_RCP85[:, :rcp_other_forcing])
run(SNEASY)
radiative_forcing_s1 = getdataframe(SNEASY, :radiativeforcing=>:rf)

# Scenario 2: update all parameters to RCP 2.6 EXCEPT leave CO2_emissions at RCP 8.5
SNEASY = MimiSNEASY.get_model()
update_param!(SNEASY, :ccm, :CO2_emissions, df_RCP85[:, :rcp_co2_emissions]) # stays at RCP 8.5
update_param!(SNEASY, :rfco2, :N₂O, df_RCP26[:, :rcp_n2o_concentration]) 
update_param!(SNEASY, :radiativeforcing, :rf_aerosol, df_RCP26[:, :rcp_aerosol_forcing])
update_param!(SNEASY, :radiativeforcing, :rf_other, df_RCP26[:, :rcp_other_forcing])
run(SNEASY)
radiative_forcing_s2 = getdataframe(SNEASY, :radiativeforcing=>:rf)

# plot scenario 1 and scenario 2 radiative forcings
years = radiative_forcing_s1[:,:time] # x-axis (years 1765-2500)
rf_scenarios = [radiative_forcing_s1[:,:rf] radiative_forcing_s2[:,:rf]]
labels = ["Scenario 1" "Scenario 2"]
rf_plot = plot(years, rf_scenarios, label = labels, legend = :topleft)
plot!(rf_plot, title = "Forcing Values for Two Scenarios", xlabel = "year", ylabel = "forcing")