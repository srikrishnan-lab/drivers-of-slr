# activate the environment
using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) # this does the same thing: Pkg.activate("/Users/ced227/Documents/Repositories/drivers-of-sea-level-rise/")
Pkg.instantiate()

# import necessary packages
using Random
using Mimi
using MimiSNEASY
using MimiBRICK
using DataFrames
using Impute
using Dates
using Distributions
using DataStructures
using XLSX
using CSVFiles # need to keep for "load" function

include("functions.jl") # include functions from other scripts

# initial set up
Random.seed!(1)
output_dir = "default"
num_samples = 1000
start_year = 1850
end_year = 2300
model_years = collect(start_year:end_year)
num_years = length(model_years)
num_params = 38

# store dataframe with info for each RCP scenario for appropriate years (contains historical data from beginning-2005)
RCP26 = scenario_df("RCP26")[in(model_years).(scenario_df("RCP26").year), :]
RCP45 = scenario_df("RCP45")[in(model_years).(scenario_df("RCP45").year), :]
RCP60 = scenario_df("RCP60")[in(model_years).(scenario_df("RCP60").year), :]
RCP85 = scenario_df("RCP85")[in(model_years).(scenario_df("RCP85").year), :]

# read in historical emissions observations (1850-2021)
historical_data = historical_emissions(start_year=start_year, end_year=end_year)

# create instance of SNEASY-BRICK
m = MimiBRICK.create_sneasy_brick(start_year=start_year, end_year=end_year) # 1 year timestep
#run(m)
#explore(m)

# generate num_samples of parameter samples for emissions curves
growth  = rand(truncated(Normal(0.004,0.0075), 0.001, Inf), num_samples)           # growth parameter
peak    = trunc.(Int64, rand(truncated(Normal(2070,25), 2030, 2200), num_samples)) # peaking time
decline = rand(truncated(Normal(0.07,0.05), 0.001, 0.2), num_samples)              # decline parameter
# add in sampling for land water storage
lw_storage =  rand(Normal(0.0003, 0.00018), num_samples, num_years)
# read in ensemble of parameter samples for SNEASY-BRICK
calibrated_params = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibrated_parameters", "parameters_subsample_sneasybrick.csv"))) # read in subsample, not full chain

# pre-allocate arrays to store results
parameters                  = zeros(Float64, num_samples, num_params)
co2_emissions               = zeros(Float64, num_samples, num_years)
radiative_forcing           = zeros(Float64, num_samples, num_years)
temperature                 = zeros(Float64, num_samples, num_years)
global_mean_sea_level_rise  = zeros(Float64, num_samples, num_years)
slr_antarctic_icesheet      = zeros(Float64, num_samples, num_years)
slr_glaciers_small_ice_caps = zeros(Float64, num_samples, num_years)
slr_greenland_icesheet      = zeros(Float64, num_samples, num_years)
slr_landwater_storage       = zeros(Float64, num_samples, num_years)
slr_thermal_expansion       = zeros(Float64, num_samples, num_years)
ocean_heat                  = zeros(Float64, num_samples, num_years)

# loop over each sample input and evaluate the model
for i = 1:num_samples

    # ----------------------------------------------- Create emissions curves with current set of samples --------------------------------------------------- #

    # isolate one set of samples and update parameters to those values
    γ_g = growth[i]
    t_peak = peak[i]
    γ_d = decline[i]

    t, gtco2 = emissions_curve(historical_data, γ_g=γ_g, t_peak=t_peak, γ_d=γ_d) # years and emissions
    gtco2 ./= 3.67 # divide by 3.67 to convert GtCO₂/yr to GtC/yr (SNEASY needs input of GtC/yr)

    # feed CO₂ emissions into SNEASY-BRICK
    update_param!(m, :ccm, :CO2_emissions, gtco2)

    # --------------------------- Match inputs for N₂O concentration, aerosol forcing, and other forcing to closest RCP scenario ---------------------------- #

    # match inputs to closest RCP scenario for each year
    emissions_df = DataFrame(Year=t, Emissions=gtco2) # need this format for match_inputs function
    n2o_concentration, aerosol_forcing, other_forcing = match_inputs(RCP26, RCP45, RCP60, RCP85, emissions_df)

    # update inputs to reflect RCP matching
    update_param!(m, :rfco2, :N₂O, n2o_concentration)
    update_param!(m, :radiativeforcing, :rf_aerosol, aerosol_forcing)
    update_param!(m, :radiativeforcing, :rf_other, other_forcing)

    # ---------------------------------------------------- Create and Update SNEASY-BRICK Parameters --------------------------------------------------------- #

    # create parameters
    lw_random_sample            = lw_storage[i,:]           # land water storage sample
    σ_temperature               = calibrated_params[i,1]    # statistical noise parameter
    σ_ocean_heat                = calibrated_params[i,2]    # statistical noise parameter
    σ_glaciers                  = calibrated_params[i,3]    # statistical noise parameter
    σ_greenland                 = calibrated_params[i,4]    # statistical noise parameter
    σ_antarctic                 = calibrated_params[i,5]    # statistical noise parameter
    σ_gmsl                      = calibrated_params[i,6]    # statistical noise parameter
    σ²_white_noise_CO₂          = calibrated_params[i,7]    # statistical noise parameter (will not use)
    ρ_temperature               = calibrated_params[i,8]    # statistical noise parameter
    ρ_ocean_heat                = calibrated_params[i,9]    # statistical noise parameter
    ρ_glaciers                  = calibrated_params[i,10]   # statistical noise parameter
    ρ_greenland                 = calibrated_params[i,11]   # statistical noise parameter
    ρ_antarctic                 = calibrated_params[i,12]   # statistical noise parameter
    ρ_gmsl                      = calibrated_params[i,13]   # statistical noise parameter
    α₀_CO₂                      = calibrated_params[i,14]   # statistical noise parameter (will not use)
    CO2_0                       = calibrated_params[i,15]
    N2O_0                       = calibrated_params[i,16]
    temperature_0               = calibrated_params[i,17]   # statistical noise parameter (will not use)
    ocean_heat_0                = calibrated_params[i,18]   # statistical noise parameter
    thermal_s0                  = calibrated_params[i,19]
    greenland_v0                = calibrated_params[i,20]
    glaciers_v0                 = calibrated_params[i,21]
    glaciers_s0                 = calibrated_params[i,22]
    antarctic_s0                = calibrated_params[i,23]
    Q10                         = calibrated_params[i,24]
    CO2_fertilization           = calibrated_params[i,25]
    CO2_diffusivity             = calibrated_params[i,26]
    heat_diffusivity            = calibrated_params[i,27]
    rf_scale_aerosol            = calibrated_params[i,28]
    climate_sensitivity         = calibrated_params[i,29]
    thermal_alpha               = calibrated_params[i,30]
    greenland_a                 = calibrated_params[i,31]
    greenland_b                 = calibrated_params[i,32]
    greenland_alpha             = calibrated_params[i,33]
    greenland_beta              = calibrated_params[i,34]
    glaciers_beta0              = calibrated_params[i,35]
    glaciers_n                  = calibrated_params[i,36]
    anto_alpha                  = calibrated_params[i,37]
    anto_beta                   = calibrated_params[i,38]
    antarctic_gamma             = calibrated_params[i,39]
    antarctic_alpha             = calibrated_params[i,40]
    antarctic_mu                = calibrated_params[i,41]
    antarctic_nu                = calibrated_params[i,42]
    antarctic_precip0           = calibrated_params[i,43]
    antarctic_kappa             = calibrated_params[i,44]
    antarctic_flow0             = calibrated_params[i,45]
    antarctic_runoff_height0    = calibrated_params[i,46]
    antarctic_c                 = calibrated_params[i,47]
    antarctic_bed_height0       = calibrated_params[i,48]
    antarctic_slope             = calibrated_params[i,49]
    antarctic_lambda            = calibrated_params[i,50]
    antarctic_temp_threshold    = calibrated_params[i,51]

    # ----- Land Water Storage ----- #
    update_param!(m, :landwater_storage, :lws_random_sample, lw_random_sample)

    # ----- Antarctic Ocean ----- #
    update_param!(m, :antarctic_ocean, :anto_α, anto_alpha)
    update_param!(m, :antarctic_ocean, :anto_β, anto_beta)

    # ----- Antarctic Ice Sheet ----- #
    update_param!(m, :antarctic_icesheet, :ais_sea_level₀, antarctic_s0)
    update_param!(m, :antarctic_icesheet, :ais_bedheight₀, antarctic_bed_height0)
    update_param!(m, :antarctic_icesheet, :ais_slope, antarctic_slope)
    update_param!(m, :antarctic_icesheet, :ais_μ, antarctic_mu)
    update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight₀, antarctic_runoff_height0)
    update_param!(m, :antarctic_icesheet, :ais_c, antarctic_c)
    update_param!(m, :antarctic_icesheet, :ais_precipitation₀, antarctic_precip0)
    update_param!(m, :antarctic_icesheet, :ais_κ, antarctic_kappa)
    update_param!(m, :antarctic_icesheet, :ais_ν, antarctic_nu)
    update_param!(m, :antarctic_icesheet, :ais_iceflow₀, antarctic_flow0)
    update_param!(m, :antarctic_icesheet, :ais_γ, antarctic_gamma)
    update_param!(m, :antarctic_icesheet, :ais_α, antarctic_alpha)
    update_param!(m, :antarctic_icesheet, :temperature_threshold, antarctic_temp_threshold)
    update_param!(m, :antarctic_icesheet, :λ, antarctic_lambda)

    # ----- Glaciers & Small Ice Caps ----- #
    update_param!(m, :glaciers_small_icecaps, :gsic_β₀, glaciers_beta0)
    update_param!(m, :glaciers_small_icecaps, :gsic_v₀, glaciers_v0)
    update_param!(m, :glaciers_small_icecaps, :gsic_s₀, glaciers_s0)
    update_param!(m, :glaciers_small_icecaps, :gsic_n, glaciers_n)

    # ----- Greenland Ice Sheet ----- #
    update_param!(m, :greenland_icesheet, :greenland_a, greenland_a)
    update_param!(m, :greenland_icesheet, :greenland_b, greenland_b)
    update_param!(m, :greenland_icesheet, :greenland_α, greenland_alpha)
    update_param!(m, :greenland_icesheet, :greenland_β, greenland_beta)
    update_param!(m, :greenland_icesheet, :greenland_v₀, greenland_v0)

    # ----- Thermal Expansion ----- #
    update_param!(m, :thermal_expansion, :te_α, thermal_alpha)
    update_param!(m, :thermal_expansion, :te_s₀, thermal_s0)

    # ----- SNEASY/DOECLIM Parameters ----- #
    update_param!(m, :doeclim, :t2co, climate_sensitivity)
    update_param!(m, :doeclim, :kappa, heat_diffusivity)
    update_param!(m, :radiativeforcing, :alpha, rf_scale_aerosol)
    update_param!(m, :model_CO₂_0, CO2_0)
    update_param!(m, :ccm, :Q10, Q10)
    update_param!(m, :ccm, :Beta, CO2_fertilization)
    update_param!(m, :ccm, :Eta, CO2_diffusivity)
    update_param!(m, :rfco2, :N₂O_0, N2O_0)

    # run SNEASY-BRICK
    run(m)

    # --------------------------------------------------- Retrieve and store output for current run -------------------------------------------------------- #

    # parameter values to be saved 
    param_vals = [γ_g, t_peak, γ_d, CO2_0, N2O_0, thermal_s0, greenland_v0, glaciers_v0, glaciers_s0, antarctic_s0, Q10, CO2_fertilization, CO2_diffusivity, heat_diffusivity, rf_scale_aerosol, climate_sensitivity, 
                 thermal_alpha, greenland_a, greenland_b, greenland_alpha, greenland_beta, glaciers_beta0, glaciers_n, anto_alpha, anto_beta, antarctic_gamma, antarctic_alpha, antarctic_mu, antarctic_nu, 
                 antarctic_precip0, antarctic_kappa, antarctic_flow0, antarctic_runoff_height0, antarctic_c, antarctic_bed_height0, antarctic_slope, antarctic_lambda, antarctic_temp_threshold]

    # write current sample to respective array
    parameters[i,:]                     = param_vals                                                # values for each parameter
    co2_emissions[i,:]                  = m[:ccm, :CO2_emissions] .* 3.67                           # total CO₂ emissions (GtCO₂/yr)
    radiative_forcing[i,:]              = m[:radiativeforcing, :rf]                                 # global radiative forcing (top of atmosphere) (W/m^2)
    temperature[i,:]                    = m[:ccm, :temp]                                            # global mean temperature anomaly (K), relative to preindustrial
    global_mean_sea_level_rise[i,:]     = m[:global_sea_level, :sea_level_rise]                     # total sea level rise from all components (m)
    slr_antarctic_icesheet[i,:]         = m[:global_sea_level, :slr_antartic_icesheet]              # sea level rise from the Antarctic ice sheet (m)
    slr_glaciers_small_ice_caps[i,:]    = m[:global_sea_level, :slr_glaciers_small_ice_caps]        # sea level rise from glaciers and small ice caps (m)
    slr_greenland_icesheet[i,:]         = m[:global_sea_level, :slr_greeland_icesheet]              # sea level rise from the Greenland ice sheet (m)
    slr_landwater_storage[i,:]          = m[:global_sea_level, :slr_landwater_storage]              # sea level rise from landwater storage (m)
    slr_thermal_expansion[i,:]          = m[:global_sea_level, :slr_thermal_expansion]              # sea level rise from thermal expansion (m)
    ocean_heat[i,:]                     = m[:doeclim, :heat_mixed] .+ m[:doeclim, :heat_interior]   # sum of ocean heat content anomaly in mixed layer and interior ocean (10²² J)

    # -------------------------------------------------- Incorporate statistical noise for current run ----------------------------------------------------- #
    
    # calculate the statistical noise
    noise_temperature   = simulate_ar1_noise(num_years, σ_temperature, ρ_temperature)
    noise_ocean_heat    = simulate_ar1_noise(num_years, σ_ocean_heat, ρ_ocean_heat)
    noise_glaciers      = simulate_ar1_noise(num_years, σ_glaciers, ρ_glaciers)
    noise_greenland     = simulate_ar1_noise(num_years, σ_greenland, ρ_greenland)
    noise_antarctic     = simulate_ar1_noise(num_years, σ_antarctic, ρ_antarctic)
    noise_gmsl          = simulate_ar1_noise(num_years, σ_gmsl, ρ_gmsl)

    # define indices for the start of statistical noise in 2010
    noise_indices = findall((in)(2010:end_year), start_year:end_year)

    # add noise to results from 2010 onwards (@view macro alters array in-place)
    @view(slr_glaciers_small_ice_caps[i,:])[noise_indices] .+= noise_glaciers[noise_indices]
    @view(slr_greenland_icesheet[i,:])[noise_indices]      .+= noise_greenland[noise_indices]
    @view(slr_antarctic_icesheet[i,:])[noise_indices]      .+= noise_antarctic[noise_indices]
    @view(global_mean_sea_level_rise[i,:])[noise_indices]  .+= noise_gmsl[noise_indices]
    @view(temperature[i,:])[noise_indices]                 .+= noise_temperature[noise_indices]
    @view(ocean_heat[i,:])[noise_indices]                  .+= noise_ocean_heat[noise_indices]

    # define baseline indices to normalize results
    temperature_norm_indices = findall((in)(1861:1880), start_year:end_year) # indices needed to normalize temperature anomalies relative to 1861-1880 mean
    sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), start_year:end_year) # indices needed to normalize sea level rise sources relative to the 1961-1990 mean
    sealevel_norm_indices_1992_2001 = findall((in)(1992:2001), start_year:end_year) # indices needed to normalize sea level rise sources relative to the 1992-2001 mean
    
    # subtract off the mean so we have results relative to a baseline
    @view(slr_glaciers_small_ice_caps[i,:]) .-= mean(slr_glaciers_small_ice_caps[i,:][sealevel_norm_indices_1961_1990])
    @view(slr_greenland_icesheet[i,:])      .-= mean(slr_greenland_icesheet[i,:][sealevel_norm_indices_1992_2001])
    @view(slr_antarctic_icesheet[i,:])      .-= mean(slr_antarctic_icesheet[i,:][sealevel_norm_indices_1992_2001])
    @view(global_mean_sea_level_rise[i,:])  .-= mean(global_mean_sea_level_rise[i,:][sealevel_norm_indices_1961_1990])
    @view(temperature[i,:])                 .-= mean(temperature[i,:][temperature_norm_indices])
    @view(ocean_heat[i,:])                  .+= ocean_heat_0 # add starting value since there is no established baseline
    
end

# save parameter names for df
param_names = ["gamma_g", "t_peak", "gamma_d", "CO2_0", "N2O_0", names(calibrated_params)[19:end]...]

# transform matrices to dataframes to improve interpretability (add parameter names/years)
parameters_df                   = DataFrame(parameters, param_names)
co2_emissions_df                = DataFrame(co2_emissions, Symbol.([model_years...]))
radiative_forcing_df            = DataFrame(radiative_forcing, Symbol.([model_years...]))
temperature_df                  = DataFrame(temperature, Symbol.([model_years...]))
global_mean_sea_level_rise_df   = DataFrame(global_mean_sea_level_rise, Symbol.([model_years...]))
slr_antarctic_icesheet_df       = DataFrame(slr_antarctic_icesheet, Symbol.([model_years...]))
slr_glaciers_small_ice_caps_df  = DataFrame(slr_glaciers_small_ice_caps, Symbol.([model_years...]))
slr_greenland_icesheet_df       = DataFrame(slr_greenland_icesheet, Symbol.([model_years...]))
slr_landwater_storage_df        = DataFrame(slr_landwater_storage, Symbol.([model_years...]))
slr_thermal_expansion_df        = DataFrame(slr_thermal_expansion, Symbol.([model_years...]))
ocean_heat_df                   = DataFrame(ocean_heat, Symbol.([model_years...]))

# export output to .csv files
save(joinpath(@__DIR__, "..", "results", "$output_dir", "parameters.csv"), parameters_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "emissions.csv"), co2_emissions_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "radiative_forcing.csv"), radiative_forcing_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "temperature.csv"), temperature_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "gmslr.csv"), global_mean_sea_level_rise_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "antarctic.csv"), slr_antarctic_icesheet_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "gsic.csv"), slr_glaciers_small_ice_caps_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "greenland.csv"), slr_greenland_icesheet_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "lw_storage.csv"), slr_landwater_storage_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "thermal_expansion.csv"), slr_thermal_expansion_df)
save(joinpath(@__DIR__, "..", "results", "$output_dir", "ocean_heat.csv"), ocean_heat_df)

#=
# create and view the netCDF file
create_netcdf(output_df)
view_nc = Dataset(joinpath(@__DIR__, "../results/slr_output.nc"), "r")
close(view_nc)
=#