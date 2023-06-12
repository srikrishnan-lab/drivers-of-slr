# this script plots the 90%, 95%, and 99% credible intervals for emissions trajectories

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots
using Random
include("../src/functions.jl")

let num_samples = 10_000
    Random.seed!(1) # set the seed

    # initialize historical observations
    historical_data = historical_emissions() # years 1850-2021

    # set up storage for years and emissions trajectories
    t = zeros(451)
    num_years = length(t)
    gtco2 = zeros(Float64, num_samples, num_years)

    # initialize plot
    p = plot(xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", legend=:topright, title="Emissions Trajectories Credible Intervals")

    # generate num_samples of parameter samples for emissions curves
    growth  = rand(truncated(Normal(0.004,0.0075), 0.001, Inf), num_samples)           # growth parameter
    peak    = trunc.(Int64, rand(truncated(Normal(2070,25), 2030, 2200), num_samples)) # peaking time
    decline = rand(truncated(Normal(0.07,0.05), 0.001, 0.2), num_samples)              # decline parameter

    # loop through samples
    for i = 1:num_samples
        # isolate one set of samples and update parameters to those values
        t[:], gtco2[i,:] = emissions_curve(historical_data, γ_g=growth[i], t_peak=peak[i], γ_d=decline[i]) # years and emissions
        #plot!(t[173:451], gtco2[i,collect(173:451)], label=:false) # add line on plot for each sample
    end

    # chop off years 1850-2021 (just get future years)
    indices = findall((in)(2022:last(t)), first(t):last(t))
    t = t[indices]
    gtco2 = gtco2[:,collect(indices)]

    # get central 90% interval for each year by picking the 5% and 95% quantiles
    gtco2_quantiles_90 = mapslices(x -> quantile(x, [0.05, 0.95]), gtco2, dims=1) # 90% interval
    gtco2_quantiles_95 = mapslices(x -> quantile(x, [0.025, 0.975]), gtco2, dims=1) # 95% interval
    gtco2_quantiles_99 = mapslices(x -> quantile(x, [0.005, 0.995]), gtco2, dims=1) # 99% interval

    # plot the credible intervals
    plot!(t, gtco2_quantiles_99[1,:], fillrange=gtco2_quantiles_99[2,:], fillalpha=0.3, alpha=0.35, color=:blue, label="99% CI")
    plot!(t, gtco2_quantiles_95[1,:], fillrange=gtco2_quantiles_95[2,:], fillalpha=0.5, alpha=0.35, color=:blue, label="95% CI")
    plot!(t, gtco2_quantiles_90[1,:], fillrange=gtco2_quantiles_90[2,:], fillalpha=0.7, alpha=0.35, color=:blue, label="90% CI")

    # add in extreme RCP scenarios for reference
    rcp26, rcp85 = rcp_emissions()
    scatter!(rcp85[:,1], rcp85[:,2], label="RCP 8.5", markersize=4, color=:black, shape=:rect) # RCP 8.5
    scatter!(rcp26[:,1], rcp26[:,2], label="RCP 2.6", markersize=5, color=:black, shape=:utriangle) # RCP 2.6

    display(p)
    #savefig(p, "/Users/ced227/Desktop/plots/emissions_CI.pdf")
    
    # calculate surprise index to check quantiles (just a sanity check)
    count = 0
    total = 0
    surprise = zeros(num_samples)
    for year in 1:length(t)
        surprise = gtco2[:,year] # num_samples length vector
        for item in surprise
            if item < gtco2_quantiles_99[1,year] # check if less than lower quantile
                count += 1
            elseif item > gtco2_quantiles_99[2,year] # check if greater than upper quantile
                count += 1
            end
            total += 1 # total count is num_samples * length(t)
        end
    end
    #println(count / total) # fraction of total count outside of specified CI
end