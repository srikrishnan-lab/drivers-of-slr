using Mimi
using MimiSNEASY
using MimiBRICK
using DataStructures

include("functions.jl")
Random.seed!(1)

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# function that establishes bandwidth for kernel density estimate (from KernelEstimator.jl)
function bwnormal(xdata::Vector)
    0.9 * min((quantile(xdata, .75) - quantile(xdata, .25)) / 1.34, std(xdata)) * length(xdata) ^ (-0.2)
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# function that samples from a Gaussian distribution near parameter values
function sample_value(;n_samples, param_vals, bw, bounds)
    # sample n_samples from a truncated Normal distribution with mean of param_vals and stdev of bw, staying inside specified bounds
    samples = rand(truncated(Normal(param_vals, bw), first(bounds), last(bounds)), n_samples)
    return samples
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

function construct_run_sneasybrick(start_year::Int, end_year::Int)

    model_years = collect(start_year:end_year)

    # Load an instance of SNEASY+BRICK model.
    # WARNING: for general use, use `m = create_sneasy_brick!(...[arguments here]...)`  instead
    m = Mimi.build(MimiBRICK.create_sneasy_brick(; rcp_scenario = "RCP85", start_year=start_year, end_year=end_year))

    # Get indices needed to normalize temperature anomalies relative to 1861-1880 mean (SNEASY+BRICK starts in 1850 by default).
    temperature_norm_indices = findall((in)(1861:1880), 1850:end_year)

    # Get indices needed to normalize all sea level rise sources.
    sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), 1850:end_year)
    sealevel_norm_indices_1992_2001 = findall((in)(1992:2001), 1850:end_year)

    # store dataframe with info for each RCP scenario for appropriate years (historical data from beginning-2005)
    RCP26 = scenario_df("RCP26")[in(model_years).(scenario_df("RCP26").year), :]
    RCP45 = scenario_df("RCP45")[in(model_years).(scenario_df("RCP45").year), :]
    RCP60 = scenario_df("RCP60")[in(model_years).(scenario_df("RCP60").year), :]
    RCP85 = scenario_df("RCP85")[in(model_years).(scenario_df("RCP85").year), :]
    # read in historical emissions observations (1850-2021)
    historical_data = historical_emissions(start_year=start_year, end_year=end_year)

    # Given user settings, create a function to run SNEASY+BRICK and return model output used for calibration.
    function run_sneasybrick!(
        param::Array{Float64,1},
        modeled_CO₂::Vector{Float64},
        modeled_oceanCO₂_flux::Vector{Float64},
        modeled_temperature::Vector{Float64},
        modeled_ocean_heat::Vector{Float64},
        modeled_glaciers::Vector{Float64},
        modeled_greenland::Vector{Float64},
        modeled_antarctic::Vector{Float64},
        modeled_thermal_expansion::Vector{Float64},
        modeled_gmsl::Vector{Float64})

        # Assign names to uncertain model and initial condition parameters for convenience.
        # Note: This assumes "param" is the full vector of uncertain parameters with the same ordering as in "create_log_posterior_sneasy_brick.jl".
        CO₂_0                    = param[16]
        N₂O_0                    = param[17]
        temperature_0            = param[18]
        ocean_heat_0             = param[19]
        thermal_s₀               = param[20]
        greenland_v₀             = param[21]
        glaciers_v₀              = param[22]
        glaciers_s₀              = param[23]
        antarctic_s₀             = param[24]
        Q10                      = param[25]
        CO₂_fertilization        = param[26]
        CO₂_diffusivity          = param[27]
        heat_diffusivity         = param[28]
        rf_scale_aerosol         = param[29]
        ECS                      = param[30]
        thermal_α                = param[31]
        greenland_a              = param[32]
        greenland_b              = param[33]
        greenland_α              = param[34]
        greenland_β              = param[35]
        glaciers_β₀              = param[36]
        glaciers_n               = param[37]
        anto_α                   = param[38]
        anto_β                   = param[39]
        antarctic_γ              = param[40]
        antarctic_α              = param[41]
        antarctic_μ              = param[42]
        antarctic_ν              = param[43]
        antarctic_precip₀        = param[44]
        antarctic_κ              = param[45]
        antarctic_flow₀          = param[46]
        antarctic_runoff_height₀ = param[47]
        antarctic_c              = param[48]
        antarctic_bedheight₀     = param[49]
        antarctic_slope          = param[50]
        antarctic_λ              = param[51]
        antarctic_temp_threshold = param[52]
        lw_random_sample         = param[53]

        #----------------------------------------------------------
        # Set SNEASY+BRICK to use sampled parameter values.
        #----------------------------------------------------------

        γ_g     = param[1]
        t_peak  = param[2]
        γ_d     = param[3]

        t, gtco2 = emissions_curve(historical_data, γ_g=γ_g, t_peak=t_peak, γ_d=γ_d, start_year=start_year, end_year=end_year) # years and emissions
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

        # ---- Land Water Storage ---- #
        update_param!(m, :landwater_storage, :lws_random_sample, fill(lw_random_sample, length(model_years)))

        # ---- Shared Parameters ---- #
        update_param!(m, :model_CO₂_0, CO₂_0) # connects to :ccm and :rfco2

        # ---- Diffusion Ocean Energy balance CLIMate model (DOECLIM) ---- #
        update_param!(m, :doeclim, :t2co,  ECS)
        update_param!(m, :doeclim, :kappa, heat_diffusivity)

        # ---- Carbon Cycle ---- #
        update_param!(m, :ccm, :Q10,     Q10)
        update_param!(m, :ccm, :Beta,    CO₂_fertilization)
        update_param!(m, :ccm, :Eta,     CO₂_diffusivity)

        # ---- Carbon Dioxide Radiative Forcing ---- #
        update_param!(m, :rfco2, :N₂O_0, N₂O_0)

        # ---- Total Radiative Forcing ---- #
        update_param!(m, :radiativeforcing, :alpha, rf_scale_aerosol)

        # ----- Antarctic Ocean ----- #
        update_param!(m, :antarctic_ocean, :anto_α, anto_α)
        update_param!(m, :antarctic_ocean, :anto_β, anto_β)

        # ----- Antarctic Ice Sheet ----- #
        update_param!(m, :antarctic_icesheet, :ais_sea_level₀,             antarctic_s₀)
        update_param!(m, :antarctic_icesheet, :ais_bedheight₀,             antarctic_bedheight₀)
        update_param!(m, :antarctic_icesheet, :ais_slope,                  antarctic_slope)
        update_param!(m, :antarctic_icesheet, :ais_μ,                      antarctic_μ)
        update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight₀, antarctic_runoff_height₀)
        update_param!(m, :antarctic_icesheet, :ais_c,                      antarctic_c)
        update_param!(m, :antarctic_icesheet, :ais_precipitation₀,         antarctic_precip₀)
        update_param!(m, :antarctic_icesheet, :ais_κ,                      antarctic_κ)
        update_param!(m, :antarctic_icesheet, :ais_ν,                      antarctic_ν)
        update_param!(m, :antarctic_icesheet, :ais_iceflow₀,               antarctic_flow₀)
        update_param!(m, :antarctic_icesheet, :ais_γ,                      antarctic_γ)
        update_param!(m, :antarctic_icesheet, :ais_α,                      antarctic_α)
        update_param!(m, :antarctic_icesheet, :temperature_threshold,      antarctic_temp_threshold)
        update_param!(m, :antarctic_icesheet, :λ,                          antarctic_λ)

        # ----- Glaciers & Small Ice Caps ----- #
        update_param!(m, :glaciers_small_icecaps, :gsic_β₀, glaciers_β₀)
        update_param!(m, :glaciers_small_icecaps, :gsic_v₀, glaciers_v₀)
        update_param!(m, :glaciers_small_icecaps, :gsic_s₀, glaciers_s₀)
        update_param!(m, :glaciers_small_icecaps, :gsic_n,  glaciers_n)

        # ----- Greenland Ice Sheet ----- #
        update_param!(m, :greenland_icesheet, :greenland_a,  greenland_a)
        update_param!(m, :greenland_icesheet, :greenland_b,  greenland_b)
        update_param!(m, :greenland_icesheet, :greenland_α,  greenland_α)
        update_param!(m, :greenland_icesheet, :greenland_β,  greenland_β)
        update_param!(m, :greenland_icesheet, :greenland_v₀, greenland_v₀)

        # ----- Thermal Expansion ----- #
        update_param!(m, :thermal_expansion, :te_α,  thermal_α)
        update_param!(m, :thermal_expansion, :te_s₀, thermal_s₀)

        # Run model.
        run(m)

        #----------------------------------------------------------
        # Calculate model output being compared to observations.
        #----------------------------------------------------------

        # Atmospheric concentration of CO₂.
        #modeled_CO₂[:] = m[:ccm, :CO₂_0]
        modeled_CO₂[:] = m[:ccm, :atmco2]

        # Global surface temperature anomaly (normalized to 1861-1880 mean with initial condition offset).
        modeled_temperature[:] = m[:doeclim, :temp] .- mean(m[:doeclim, :temp][temperature_norm_indices]) .+ temperature_0

        # Ocean carbon flux (Note: timesteps cause last `atm_oc_flux` value to equal `missing`, so exclude it here).
        modeled_oceanCO₂_flux[1:end-1] = m[:ccm, :atm_oc_flux][1:end-1]

        # Ocean heat content (with initial condition offset).
        modeled_ocean_heat[:] = m[:doeclim, :heat_mixed] .+ m[:doeclim, :heat_interior] .+ ocean_heat_0

        # Glaciers and small ice caps (normalized relative to 1961-1990 mean).
        modeled_glaciers[:] = m[:glaciers_small_icecaps, :gsic_sea_level] .- mean(m[:glaciers_small_icecaps, :gsic_sea_level][sealevel_norm_indices_1961_1990])

        # Greenland ice sheet (normalized relative to 1992-2001 ten year period to work with pooled data that includes IMBIE observations).
        modeled_greenland[:] = m[:greenland_icesheet, :greenland_sea_level] .- mean(m[:greenland_icesheet, :greenland_sea_level][sealevel_norm_indices_1961_1990])

        # Antarctic ice sheet (normalized relative to 1992-2001 ten year period to work with IMBIE data).
        modeled_antarctic[:] = m[:antarctic_icesheet, :ais_sea_level] .- mean(m[:antarctic_icesheet, :ais_sea_level][sealevel_norm_indices_1992_2001])

        # Sea level contribution from thermal expansion (calibrating to observed trends, so do not need to normalize).
        modeled_thermal_expansion[:] = m[:thermal_expansion, :te_sea_level]

		# Global mean sea level rise (normalize realtive to 1961-1990 mean).
		modeled_gmsl[:] = m[:global_sea_level, :sea_level_rise] .- mean(m[:global_sea_level, :sea_level_rise][sealevel_norm_indices_1961_1990])

        # Return results.
        return
    end

    # Return run model function.
    return run_sneasybrick!
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

function model_ensemble(M::Matrix{Float64}; start_year=1850, end_year=2300)
    # initial set up
    num_samples = size(M,2)
    model_years = collect(start_year:end_year)
    num_years = length(model_years)

    # pre-allocate vectors to store model output
    modeled_CO₂               = zeros(num_years)
    modeled_oceanCO₂_flux     = zeros(num_years)
    modeled_temperature       = zeros(num_years)
    modeled_ocean_heat        = zeros(num_years)
    modeled_glaciers          = zeros(num_years)
    modeled_greenland         = zeros(num_years)
    modeled_antarctic         = zeros(num_years)
    modeled_thermal_expansion = zeros(num_years)
    modeled_gmsl              = zeros(num_years)
    gmsl_out                  = zeros(num_samples)

    # create function to run SNEASY-BRICK
    run_sneasy_brick! = construct_run_sneasybrick(start_year, end_year) 

    # subset df to get values for desired year
    yr_index = findfirst(end_year .∈ model_years) # gets index for column of end_year

    # loop over each sample input and evaluate the model
    for i = 1:num_samples
        # ----------------------------------------------- Create emissions curves with current set of samples --------------------------------------------------- #
        
        # isolate one set of samples and update parameters to those values
        p = M[:, i]

        # ---------------------------------------------------- Create and Update SNEASY-BRICK Parameters --------------------------------------------------------- #

        # create parameters
        σ_temperature               = p[4]  # statistical noise parameter
        σ_ocean_heat                = p[5]  # statistical noise parameter
        σ_glaciers                  = p[6]  # statistical noise parameter
        σ_greenland                 = p[7]  # statistical noise parameter
        σ_antarctic                 = p[8]  # statistical noise parameter
        σ_gmsl                      = p[9]  # statistical noise parameter
        ρ_temperature               = p[10] # statistical noise parameter
        ρ_ocean_heat                = p[11] # statistical noise parameter
        ρ_glaciers                  = p[12] # statistical noise parameter
        ρ_greenland                 = p[13] # statistical noise parameter
        ρ_antarctic                 = p[14] # statistical noise parameter
        ρ_gmsl                      = p[15] # statistical noise parameter
   
        run_sneasy_brick!(p, modeled_CO₂, modeled_oceanCO₂_flux, modeled_temperature, modeled_ocean_heat,
        modeled_glaciers, modeled_greenland, modeled_antarctic, modeled_thermal_expansion, modeled_gmsl)

        # -------------------------------------------------- Incorporate statistical noise for current run ----------------------------------------------------- #
        
        # calculate the statistical noise
        noise_gmsl = simulate_ar1_noise(num_years, σ_gmsl, ρ_gmsl)

        # define indices for the start of statistical noise in 2010
        noise_indices = findall((in)(2010:end_year), start_year:end_year)

        # add noise to results from 2010 onwards 
        simulated_gmsl = modeled_gmsl
        simulated_gmsl[noise_indices] += noise_gmsl[noise_indices]

        # define baseline indices to normalize results relative to the 1961-1990 mean
        sealevel_norm_indices_1961_1990 = findall((in)(1961:1990), start_year:end_year)

        # subtract off the mean so we have results relative to the baseline
        simulated_gmsl .-= mean(simulated_gmsl[sealevel_norm_indices_1961_1990])

        # write current sample to array
        gmsl_out[i] = simulated_gmsl[yr_index]   # total sea level rise from all components (m)
    end

    return gmsl_out
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# function that takes in an input matrix M where each row is a set of parameters, and then returns a vector of GMSLR for a given year
function brick_run(M::Matrix{Float64}; yr=2100)
    # intialize values
    start_year = 1850

    # function returns a df of global mean sea level rise values
    gmslr = model_ensemble(M, start_year=start_year, end_year=yr)

    return (gmslr .- mean(gmslr)) / std(gmslr) # return centered GMSLR for the specified year (vector with length n_samples)
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #