# this script is used to visualize low, medium, & high emissions scenarios
# tested different growth, t_peak, & decline parameters to see best fit

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots
include("../src/functions.jl")

# initialize values
start_year = 1850
end_year = 2300
model_years = collect(start_year:end_year)
num_years = length(model_years)

# initialize historical observations
historical_data = historical_emissions() # years 1850-2021

n = plot(title="Selected CO₂ Emissions Trajectories", xlabel="Year", ylabel="Total CO₂ Emissions (GtCO₂/yr)", legend=:topleft)

# low emissions
t_low, gtco2_low = emissions_curve(historical_data; γ_g=0.006, γ_d=0.15, t_peak=2030) # fast decline and early peak
plot!(t_low, gtco2_low, label="Low")
# med emissions
t_med, gtco2_med = emissions_curve(historical_data; γ_g=0.009, γ_d=0.05, t_peak=2075) # middle of the road
plot!(t_med, gtco2_med, label="Medium")
# high emissions
t_high, gtco2_high = emissions_curve(historical_data; γ_g=0.012, γ_d=0.02, t_peak=2120) # slow decline and late peak
plot!(t_high, gtco2_high, label="High")

# rcp scenarios
rcp26, rcp85 = rcp_emissions()
scatter!(n, rcp26[:,1], rcp26[:,2], label="RCP 2.6", color=:black, markersize=3, shape=:rect)
scatter!(n, rcp85[:,1], rcp85[:,2], label="RCP 8.5", color=:black, markersize=4, shape=:utriangle)

display(n)