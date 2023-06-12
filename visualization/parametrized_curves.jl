# this script is used to plot 100 samples of emissions trajectories to visualize behavior

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots
include("../src/functions.jl")

# ------------------------------------------- Plot 100 samples of emissions trajectories --------------------------------------------- #
let num_samples = 100
    Random.seed!(1) # set the seed

    # initialize historical observations
    historical_data = historical_emissions() # years 1850-2021

    # data for extreme RCP scenarios
    rcp26, rcp85 = rcp_emissions() # units in GtCO₂

    # initiate a plot
    p = scatter(xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", title="Selection of 100 Emissions Trajectories", legend=:topright)
                #tickfontsize=20, labelfontsize=22, legendfontsize=12, titlefontsize=26, size=(1000,700), margin=5Plots.mm)
    
    # loop through samples
    for i in 1:num_samples
        # sample parameters for current iteration
        γ_g     = rand(truncated(Normal(0.004,0.0075), 0.001, Inf))           # growth parameter
        t_peak  = trunc.(Int64, rand(truncated(Normal(2070,25), 2030, 2200))) # peaking time
        γ_d     = rand(truncated(Normal(0.07,0.05), 0.001, 0.2))              # decline parameter
        
        # years and emissions for one iteration
        t, gtco2 = emissions_curve(historical_data, γ_g=γ_g, t_peak=t_peak, γ_d=γ_d)
            
        # add the curve for this iteration to the plot from the year 2022 onward
        idx = findfirst(2022 .∈ first(t):last(t))
        plot!(t[idx:end], gtco2[idx:end], label=:false)
    end

    # historical data & points that we want to include in emissions trajectories
    indices = findall((in)(2000:last(historical_data.year)), first(historical_data.year):last(historical_data.year))
    scatter!(historical_data.year[indices], historical_data.rcp_co2_emissions[indices], label="Historical Data", markersize=3, color=:black) # observed data
    scatter!(rcp85[:,1], rcp85[:,2], label="RCP 8.5", markersize=4, color=:black, shape=:rect) # RCP 8.5
    scatter!(rcp26[:,1], rcp26[:,2], label="RCP 2.6", markersize=5, color=:black, shape=:utriangle) # RCP 2.6

    # display and save the final plot
    display(p)
    #savefig(p, "/Users/ced227/Desktop/plots/100_curves.pdf")
end