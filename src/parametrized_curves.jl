using Plots
using Polynomials
using Distributions

include("exploratory_plots.jl")

# ------------------------------------------- Create baseline curve for the a likely emissions trajectory --------------------------------------------- #

emissions = [[2030,41],
            [2032,42.5],
            [2035,43],
            [2050,35],
            [2100,17],
            [2150,11],
            [2200,9],
            [2250,8]]

# plot and generate best-fit equations
plot(first.(emissions), last.(emissions), xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", title="Base Curve", legend=false)
growth_eqn = fit(first.(emissions)[1:3], last.(emissions)[1:3], 2) 
decline_eqn = fit(first.(emissions)[3:8], last.(emissions)[3:8], 2)
scatter!(growth_eqn, 2030, 2035)
scatter!(decline_eqn, 2035, 2250)

# ------------------------------------------- Develop piecewise function for parametrized curves --------------------------------------------- #

let t = zeros(Int64, 280), gtco2 = zeros(280) # allocate space for years and emissions (range from years 2021-2300)

    # 36.3 GtCO₂ emitted in 2021 (start year)
    t[1] = 2021
    gtco2[1] = 36.3
    gtco2_tpeak = 0 # need this in scope for later

    # points that we want to include in emissions trajectories
    p = scatter(first.(emissions), last.(emissions), xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", title="Emissions Trajectories", label="Baseline Data", legend=:topleft)
    scatter!(rcp26[:,1], rcp26[:,2], label="RCP 2.6") # RCP 2.6
    scatter!(rcp85[:,1], rcp85[:,2], label="RCP 8.5") # RCP 8.5

    # sample parameters for each run
    for run in 1:100
        γ_g    = rand(Uniform(0.001,0.012)) # growth parameter  (original = 0.0089), rand(Uniform(0,0.02))
        γ_d    = rand(Uniform(0.005,0.13)) # decline parameter, rand(Uniform(0.005,0.126))
        t_peak = trunc(Int64, rand(Uniform(2030,2200))) # peaking time     (original = 2115),   rand(Uniform(2030,2200))

        for i in 2:length(t)
            t[i] = t[i-1] + 1 # fill in current year to the t array
            # before peaking time: quadratic increase
            if t[i] <= t_peak
                gtco2[i] = gtco2[i-1] + γ_g * (t_peak - t[i]) # calculate emissions for current year
                if t[i] == t_peak
                    gtco2_tpeak = gtco2[i-1] + γ_g * (t_peak - t[i]) # save value of emissions at t_peak for scaling in logistic function
                end
            # after peaking time: logistic decrease
            elseif t[i] > t_peak
                gtco2[i] = gtco2[i-1] - ((2 * gtco2_tpeak * γ_d * exp(γ_d * (t[i] - t_peak))) / (exp(γ_d * (t[i] - t_peak)) + 1)^2) # calculate emissions for current year
            end
        end
        
    # add the curve for this run to the plot
    plot!(t, gtco2, label=false)
      
    end

    # display the final plot
    display(p)

end

# plots growth_eqn and decline_eqn fluidly
plot(x -> piecewise(x, t_peak=5), -9, 220, xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", legend=:bottomleft)

# ------------------------------------------- Add in extreme RCP Scenarios to ensure they are covered in generated curves --------------------------------------------- #

# anthropogenic total CO2 emissions (PgC/yr)
rcp26 = [2000 8.03;
        2010 9.70;
        2020 9.97;
        2030 8.00;
        2040 5.30;
        2050 3.50;
        2060 2.10;
        2070 0.81;
        2080 0.16;
        2090 -0.23;
        2100 -0.42]

# anthropogenic total CO2 emissions (PgC/yr)
rcp85 = [2000 8.03;
        2010 9.98;
        2020 12.28;
        2030 14.53;
        2040 17.33;
        2050 20.61;
        2060 23.83;
        2070 26.17;
        2080 27.60;
        2090 28.44;
        2100 28.77]

# multiply by 3.67 to convert PgC/yr to GtCO₂/yr (1 PgC = 1 GtC)
rcp26[:,2] = rcp26[:,2] * 3.67
rcp85[:,2] = rcp85[:,2] * 3.67

# add RCP scenarios to plot as lower/upper bounds
scatter(rcp26[:,1], rcp26[:,2], label="RCP 2.6", title="RCP 2.6")
scatter!(rcp85[:,1], rcp85[:,2], label="RCP 8.5", title="RCP 8.5", legend=:topleft)

# fit equations and plot (don't need to fit equations)
rcp26_eqn = fit(rcp26[:,1], rcp26[:,2], 3)
plot!(rcp26_eqn, extrema(rcp26[:,1])..., label="best fit 2.6")

rcp85_eqn = fit(rcp85[:,1], rcp85[:,2], 3)
plot!(rcp85_eqn, extrema(rcp85[:,1])..., label="best fit 8.5")