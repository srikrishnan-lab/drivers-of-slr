# activate the environment
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots
using Polynomials
using Distributions

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
plot(first.(emissions), last.(emissions), xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", title="Base Curve", legend=false)
growth_eqn = Polynomials.fit(first.(emissions)[1:3], last.(emissions)[1:3], 2) 
decline_eqn = Polynomials.fit(first.(emissions)[3:8], last.(emissions)[3:8], 2)
scatter!(growth_eqn, 2030, 2035)
scatter!(decline_eqn, 2035, 2250)

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

# ------------------------------------------- Develop piecewise function for parametrized curves --------------------------------------------- #

let t = zeros(Int64, 280), gtco2 = zeros(280), gtco2_tpeak = 0 # allocate space for years and emissions (range from years 2021-2300)

    # 36.3 GtCO₂ emitted in 2021 (start year)
    t[1] = 2021
    gtco2[1] = 36.3

    # points that we want to include in emissions trajectories
    p = scatter(first.(emissions), last.(emissions), xlabel="Year", ylabel="CO₂ Emissions (GtCO₂)", title="Emissions Trajectories", label="Baseline Data", legend=:topleft)
    scatter!(rcp26[:,1], rcp26[:,2], label="RCP 2.6") # RCP 2.6
    scatter!(rcp85[:,1], rcp85[:,2], label="RCP 8.5") # RCP 8.5

    for run in 1:100 # loop through runs

        # sample parameters for each run
        γ_g    = rand(Uniform(0.001,0.012))             # growth parameter      (original = 0.0089), rand(Uniform(0,0.02))
        γ_d    = rand(Uniform(0.005,0.13))              # decline parameter     rand(Uniform(0.005,0.126))
        t_peak = trunc(Int64, rand(Uniform(2030,2200))) # peaking time          (original = 2115), rand(Uniform(2030,2200))

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
    plot!(t, gtco2, label=false)
      
    end

    # display the final plot
    display(p)

end