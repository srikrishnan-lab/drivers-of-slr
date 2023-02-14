using CSVFiles
using DataFrames
using Mimi
using Plots
using MimiSNEASY

include("functions.jl")

#---------------------------------------------------------------------- Plot 1 ----------------------------------------------------------------------#

# create dataframe for each RCP scenario
df_RCP26 = scenario_df("RCP26")
df_RCP45 = scenario_df("RCP45")
df_RCP60 = scenario_df("RCP60")
df_RCP85 = scenario_df("RCP85")

# years and labels for plots
years = [1765:1:2500;]
labels = ["RCP 2.6" "RCP 4.5" "RCP 6.0" "RCP 8.5"]

# plot rcp_co2_emissions
all_rcp_co2_emissions = [df_RCP26[:,:rcp_co2_emissions] df_RCP45[:,:rcp_co2_emissions] df_RCP60[:,:rcp_co2_emissions] df_RCP85[:,:rcp_co2_emissions]]
emissions_plot = plot(years, all_rcp_co2_emissions, label = labels, title = "rcp_co2_emissions", xlabel = "Year", ylabel = "CO₂ Emissions (GtCO₂/yr)")

# plot rcp_n2o_concentration
all_rcp_n2o_concentration = [df_RCP26[:,:rcp_n2o_concentration] df_RCP45[:,:rcp_n2o_concentration] df_RCP60[:,:rcp_n2o_concentration] df_RCP85[:,:rcp_n2o_concentration]]
n2o_plot = plot(years, all_rcp_n2o_concentration, label = labels, legend = :topleft, title = "rcp_n2o_concentration", xlabel = "Year", ylabel = "N₂O Concentration (ppb)")

# plot rcp_aerosol_forcing
all_rcp_aerosol_forcing = [df_RCP26[:,:rcp_aerosol_forcing] df_RCP45[:,:rcp_aerosol_forcing] df_RCP60[:,:rcp_aerosol_forcing] df_RCP85[:,:rcp_aerosol_forcing]]
aerosol_plot = plot(years, all_rcp_aerosol_forcing, label = labels, title = "rcp_aerosol_forcing", xlabel = "Year", ylabel = "Aerosol Forcing (W/m₂)")

# plot rcp_other_forcing
all_rcp_other_forcing = [df_RCP26[:,:rcp_other_forcing] df_RCP45[:,:rcp_other_forcing] df_RCP60[:,:rcp_other_forcing] df_RCP85[:,:rcp_other_forcing]]
other_forcing_plot = plot(years, all_rcp_other_forcing, label = labels, legend = :topleft, title = "rcp_other_forcing", xlabel = "Year", ylabel = "Other Forcing (W/m₂)")

# display all plots
all_plots = plot(emissions_plot, n2o_plot, aerosol_plot, other_forcing_plot, size=(1300,1100))
display(all_plots)

#---------------------------------------------------------------------- Plot 2 ----------------------------------------------------------------------#

# Scenario 1: update all parameters to RCP 8.5
SNEASY = MimiSNEASY.get_model()
update_param!(SNEASY, :ccm, :CO2_emissions, df_RCP85[:, :rcp_co2_emissions]./3.67) # SNEASY CO₂ emissions are in GtC/yr, so need to convert GtCO₂/yr -> GtC/yr to match units (divide by 3.67)
update_param!(SNEASY, :rfco2, :N₂O, df_RCP85[:, :rcp_n2o_concentration]) 
update_param!(SNEASY, :radiativeforcing, :rf_aerosol, df_RCP85[:, :rcp_aerosol_forcing])
update_param!(SNEASY, :radiativeforcing, :rf_other, df_RCP85[:, :rcp_other_forcing])
run(SNEASY)
radiative_forcing_s1 = getdataframe(SNEASY, :radiativeforcing=>:rf)

# Scenario 2: update all parameters to RCP 2.6 EXCEPT leave CO2_emissions at RCP 8.5
SNEASY = MimiSNEASY.get_model()
update_param!(SNEASY, :ccm, :CO2_emissions, df_RCP85[:, :rcp_co2_emissions]./3.67) # stays at RCP 8.5
update_param!(SNEASY, :rfco2, :N₂O, df_RCP26[:, :rcp_n2o_concentration]) 
update_param!(SNEASY, :radiativeforcing, :rf_aerosol, df_RCP26[:, :rcp_aerosol_forcing])
update_param!(SNEASY, :radiativeforcing, :rf_other, df_RCP26[:, :rcp_other_forcing])
run(SNEASY)
radiative_forcing_s2 = getdataframe(SNEASY, :radiativeforcing=>:rf)

# plot Scenario 1 and Scenario 2 radiative forcings
years = radiative_forcing_s1[:,:time] # x-axis (years 1765-2500)
rf_scenarios = [radiative_forcing_s1[:,:rf] radiative_forcing_s2[:,:rf]]
labels = ["Scenario 1" "Scenario 2"]
rf_plot = plot(years, rf_scenarios, label = labels, legend = :topleft, title = "Forcing Values for Two Scenarios", xlabel = "Year", ylabel = "Radiative Forcing (W/m₂)")
display(rf_plot)