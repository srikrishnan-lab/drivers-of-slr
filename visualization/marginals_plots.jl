# this script plots the priors and marginal posteriors for the Antarctic and Greenland parameters in BRICK.
# Code for priors was taken from MimiBRICK.jl/src/calibration/create_log_posteriors/create_log_posterior_sneasybrick.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using NetCDF
using KernelDensity
using Distributions
using StatsPlots
using CSVFiles
using DataFrames
using Measures

# --------------------------------------- #
# --------------- Priors ---------------- #
# --------------------------------------- #

# this function calculates kernel density estimates with truncated bounds
function truncated_kernel(data, lower_bound, upper_bound)
    # Estimate kernel density without specifying boundaries to avoid wrap-around at margins.
    kde_0 = kde(data)
    # Get number of values in KDE grid.
    n = length(kde_0.x)
    # Construct range of values that spans the KDE.
    kde_range = LinRange(minimum(kde_0.x), maximum(kde_0.x), n)
    # Get indices correponding to kde_range that represent desired truncated limits.
    lower_index = searchsortedfirst(kde_range, lower_bound)
    upper_index = searchsortedlast(kde_range, upper_bound)
    # Construct a truncated KDE object.
    truncated_kde = UnivariateKDE(kde_0.x[lower_index:upper_index], kde_0.density[lower_index:upper_index])
    # Interpolate the truncated KDE object for efficiency.
    # Note: "If you are making multiple calls to pdf, it will be more efficient to construct an intermediate InterpKDE to store the interpolation structure:"
    truncated_kde_interp = InterpKDE(truncated_kde)
    # Return truncated and interpolated KDE.
    return truncated_kde_interp
end

# Load required data to create Antarctic ice sheet informative priors (posterior parameters from previous calibration to paleo data).
antarctic_paleo_file   = joinpath(@__DIR__, "..", "data", "DAISfastdyn_calibratedParameters_gamma_29Jan2017.nc")
antarctic_paleo_params = convert(Array{Float64,2}, ncread(antarctic_paleo_file, "DAIS_parameters"))'[:,1:15]

# Calculate upper and lower bounds for Antarctic ice sheet parameters (min/max values from paleo calibration).
antarctic_lower_bound = vec(minimum(antarctic_paleo_params, dims=1))
antarctic_upper_bound = vec(maximum(antarctic_paleo_params, dims=1))

# get Antarctic priors
prior_anto_α         = truncated_kernel(antarctic_paleo_params[:,1],  antarctic_lower_bound[1],  antarctic_upper_bound[1])
prior_anto_β         = truncated_kernel(antarctic_paleo_params[:,2],  antarctic_lower_bound[2],  antarctic_upper_bound[2])
prior_γ              = truncated_kernel(antarctic_paleo_params[:,3],  antarctic_lower_bound[3],  antarctic_upper_bound[3])
prior_α              = truncated_kernel(antarctic_paleo_params[:,4],  antarctic_lower_bound[4],  antarctic_upper_bound[4])
prior_μ              = truncated_kernel(antarctic_paleo_params[:,5],  antarctic_lower_bound[5],  antarctic_upper_bound[5])
prior_ν              = truncated_kernel(antarctic_paleo_params[:,6],  antarctic_lower_bound[6],  antarctic_upper_bound[6])
prior_precip₀        = truncated_kernel(antarctic_paleo_params[:,7],  antarctic_lower_bound[7],  antarctic_upper_bound[7])
prior_κ              = truncated_kernel(antarctic_paleo_params[:,8],  antarctic_lower_bound[8],  antarctic_upper_bound[8])
prior_flow₀          = truncated_kernel(antarctic_paleo_params[:,9],  antarctic_lower_bound[9],  antarctic_upper_bound[9])
prior_runoff_height₀ = truncated_kernel(antarctic_paleo_params[:,10], antarctic_lower_bound[10], antarctic_upper_bound[10])
prior_c              = truncated_kernel(antarctic_paleo_params[:,11], antarctic_lower_bound[11], antarctic_upper_bound[11])
prior_bedheight₀     = truncated_kernel(antarctic_paleo_params[:,12], antarctic_lower_bound[12], antarctic_upper_bound[12])
prior_slope          = truncated_kernel(antarctic_paleo_params[:,13], antarctic_lower_bound[13], antarctic_upper_bound[13])
prior_λ              = truncated_kernel(antarctic_paleo_params[:,14], antarctic_lower_bound[14], antarctic_upper_bound[14])
prior_temp_threshold = truncated_kernel(antarctic_paleo_params[:,15], antarctic_lower_bound[15], antarctic_upper_bound[15])
prior_antarctic_s₀   = Uniform(-0.04755, 0.05585)

# get Greenland priors
prior_greenland_v₀   = Uniform(7.16, 7.56)
prior_greenland_a    = Uniform(-4.0, -0.001)
prior_greenland_b    = Uniform(5.888, 8.832)
prior_greenland_α    = Uniform(0.0, 0.001)
prior_greenland_β    = Uniform(0.0, 0.001)

priors_1   = [prior_anto_α, prior_anto_β, prior_γ, prior_α, prior_μ, prior_ν, prior_precip₀, prior_κ, prior_flow₀,
              prior_runoff_height₀, prior_c, prior_bedheight₀, prior_slope, prior_λ, prior_temp_threshold]
priors_2   = [prior_antarctic_s₀, prior_greenland_v₀, prior_greenland_a, 
              prior_greenland_b, prior_greenland_α, prior_greenland_β]
all_priors = [priors_1..., priors_2...]

# --------------------------------------- #
# -------- Marginal Posteriors ---------- #
# --------------------------------------- #

# read in subsample of MCMC parameters (10,000 samples)
mcmc_params = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibrated_parameters", "parameters_subsample_sneasybrick.csv")))

all_posteriors = ["anto_alpha", "anto_beta", "antarctic_gamma", "antarctic_alpha", "antarctic_mu", "antarctic_nu", "antarctic_precip0",
                  "antarctic_kappa", "antarctic_flow0", "antarctic_runoff_height0", "antarctic_c", "antarctic_bed_height0", "antarctic_slope",
                  "antarctic_lambda", "antarctic_temp_threshold", "antarctic_s0", "greenland_v0", "greenland_a", "greenland_b",
                  "greenland_alpha", "greenland_beta"]

# --------------------------------------- #
# -- Plot Priors & Marginal Posteriors -- #
# --------------------------------------- #

let
    all_plots = []
    i = 1
    j = 1
    for (prior, posterior) in zip(all_priors, all_posteriors)
        # intialize plot
        p1 = plot(title="$posterior", xlabel="Parameter Value", ylabel="Density")
        # plot marginal posterior
        density!(p1, mcmc_params[:,posterior], label="Posterior", linewidth=2.5)
        # create a tuple with min/max values for the parameter
        bounds = extrema(mcmc_params[:, posterior])
        # create values for x-axis
        s = 0.0001 # step between lower and upper value
        x = first(bounds):s:last(bounds)
        # add priors to plot (change plot function call depending on prior type)
        if prior ∈ priors_1 # if the prior is an InterpKDE type
            xs = [x, 0:s:2.5, x, # tweak bounds if needed
                  x, 6:s:15, x,
                  x, x, x,
                  x, x, 730:s:840,
                  x, x, x]
            density!(p1, pdf(prior, xs[i]), xs[i], label="Prior", linewidth=2.5, linestyle=:dash)
            i += 1 # increment iterator
        elseif prior ∈ priors_2 # if the prior is a uniform distribution
            xs = [-0.06:s:0.1, 7.1:s:7.74, -4.5:s:0.5,
                5.7:s:9.5, -0.0005:s:0.0015, -0.00002:s:0.0011]
            plot!(p1, xs[j], pdf.(prior, xs[j]), trim=false, label="Prior", linewidth=2.5, linestyle=:dash)
            j += 1 # increment iterator
        end
        push!(all_plots, p1) # add current plot to vector
    end
    p2 = plot(all_plots..., layout=(7,3), size=(1800,2500), titlefontsize=17, tickfontsize=11, labelfontsize=11, plot_titlefontsize=30,
            left_margin=15mm, bottom_margin=5mm, plot_title="Priors and Marginal Posteriors for Selected BRICK Parameters")
    display(p2)
    #savefig(p2, "/Users/ced227/Desktop/plots/marginals.pdf")
end