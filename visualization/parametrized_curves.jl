# this script is used to plot a sample of 100 emissions trajectories to visualize behavior

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots
using Polynomials
using Random
include("../src/functions.jl")

# ------------------------------------------- Create baseline curve for a likely emissions trajectory --------------------------------------------- #

# estimated from Rennert et al., 2022 - Figure 1(c)
emissions = [[2030,41],
            [2032,42.5],
            [2035,43],
            [2050,35],
            [2100,17],
            [2150,11],
            [2200,9],
            [2250,8]] # units in GtCO₂

# plot and generate best-fit equations
plot(first.(emissions), last.(emissions), xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", title="Base Curve", legend=:false)
growth_eqn = Polynomials.fit(first.(emissions)[1:3], last.(emissions)[1:3], 2) 
decline_eqn = Polynomials.fit(first.(emissions)[3:8], last.(emissions)[3:8], 2)
scatter!(growth_eqn, 2030, 2035)
scatter!(decline_eqn, 2035, 2250)

# ------------------------------------------- Develop piecewise function for parametrized curves --------------------------------------------- #

let t = zeros(Int64, 280), gtco2 = zeros(280), gtco2_tpeak = 0 # allocate space for years and emissions (range from years 2021-2300)

    Random.seed!(1) # set the seed

    # initialize historical observations
    historical_data = historical_emissions() # years 1850-2021

    # data for extreme RCP scenarios
    rcp26, rcp85 = rcp_emissions() # units in GtCO₂

    # 41.13 GtCO₂ emitted in 2021 (start year)
    t[1] = last(historical_data.year)
    gtco2[1] = last(historical_data.rcp_co2_emissions)

    # initiate a plot
    p = scatter(xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", title="Sample of 100 Emissions Trajectories", legend=:topright,
                tickfontsize=20, labelfontsize=22, legendfontsize=12, titlefontsize=26, size=(1000,700), margin=5Plots.mm)

    for run in 1:100 # loop through runs

        # sample parameters for each run
        γ_g     = rand(truncated(Normal(0.004,0.0075), 0.001, Inf))           # growth parameter
        t_peak  = trunc.(Int64, rand(truncated(Normal(2070,25), 2030, 2200))) # peaking time
        γ_d     = rand(truncated(Normal(0.07,0.05), 0.001, 0.2))              # decline parameter

        for i in 2:length(t) # loop through years
            
            t[i] = t[i-1] + 1 # fill in current year to the t array

            # before peaking time: quadratic increase
            if t[i] <= t_peak
                gtco2[i] = gtco2[i-1] + γ_g * (t_peak - t[i]) # calculate emissions for current year
                
                # at peaking time: save value
                if t[i] == t_peak
                    gtco2_tpeak = gtco2[i-1] + γ_g * (t_peak - t[i]) # save value of emissions at t_peak for scaling in logistic function
                end

            # after peaking time: logistic decrease
            elseif t[i] > t_peak
                gtco2[i] = gtco2[i-1] - ((2 * gtco2_tpeak * γ_d * exp(γ_d * (t[i] - t_peak))) / (exp(γ_d * (t[i] - t_peak)) + 1)^2) # calculate emissions for current year
            end

        end
        
    # add the curve for this run to the plot
    plot!(t, gtco2, label=:false)
      
    end

    # historical data & points that we want to include in emissions trajectories
    idx = findall((in)(2000:last(historical_data.year)), first(historical_data.year):last(historical_data.year))
    scatter!(historical_data.year[idx], historical_data.rcp_co2_emissions[idx], label="Historical Data", markersize=4, color=:black) # observed data
    scatter!(rcp85[:,1], rcp85[:,2], label="RCP 8.5", markersize=5, color=:black, shape=:rect) # RCP 8.5
    scatter!(rcp26[:,1], rcp26[:,2], label="RCP 2.6", markersize=6, color=:black, shape=:utriangle) # RCP 2.6

    # display/save the final plot
    display(p)
    #savefig(p, "/Users/ced227/Desktop/plots/output.png")

end