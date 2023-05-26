using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using QuasiMonteCarlo
include("gsa_functions.jl")

Random.seed!(1)
n_samples = 10_000 # define size of ensemble
output_dir = "default"

#----------------------------------------------- Create design matrices A1 and B1 (for emissions parameters) ------------------------------------------------------------------------------------#

# set bounds for the 3 emissions parameters (growth, t_peak, decline)
lb = zeros(3) # lower bound (0)
ub = ones(3) # upper bound (1)

# create design matrices for emissions parameters (A1 and B1)
A1,B1 = QuasiMonteCarlo.generate_design_matrices(n_samples, lb, ub, SobolSample()) # all values in range 0 to 1

# define distributions for emissions parameters
γ_g_dist     = truncated(Normal(0.004,0.0075), 0.001, Inf)  # growth parameter
t_peak_dist  = truncated(Normal(2070,25), 2030, 2200)       # peaking time
γ_d_dist     = truncated(Normal(0.07,0.05), 0.001, 0.2)     # decline parameter

# stratefied sensitivity analysis runs
#t_peak_dist  = truncated(Normal(2070,25), 2030, 2050)       # early peaking time (before 2050)
#t_peak_dist  = truncated(Normal(2070,25), 2050, 2100)       # middle peaking time (2050-2100)
#t_peak_dist  = truncated(Normal(2070,25), 2100, 2200)       # late peaking time (after 2100)

# combine emissions parmaters' distributions into a matrix
emissions_dist = hcat(γ_g_dist, t_peak_dist, γ_d_dist)

# convert the 0-1 quantiles into the actual values for growth, peaking, and decline
for i in 1:3 # loop through the three emissions parameters
    A1[i,:] = quantile(emissions_dist[i], A1[i,:])
    B1[i,:] = quantile(emissions_dist[i], B1[i,:])
end

# truncate the peaking years from Float to Integer
A1[2,:] = trunc.(Int64, A1[2,:])
B1[2,:] = trunc.(Int64, B1[2,:])

#----------------------------------------------- Create design matrices A2 and B2 (for Earth system parameters) ----------------------------------------------------------------------#

# read in subsample of MCMC parameters (10,000 samples)
mcmc_params = DataFrame(load(joinpath(@__DIR__, "..", "data", "calibrated_parameters", "parameters_subsample_sneasybrick.csv"))) # read in subsample

# get max and min values for each posterior distribution's samples
bounds = [mapslices(extrema, Matrix(mcmc_params), dims=1)...] # produces a vector of tuples with min/max values for each parameter

# fit KDEs to each marginal posterior (we need output of bandwidth)
bandwidth = mapslices(bwnormal, Matrix(mcmc_params); dims=1)

# initialize storage for A2 and B2 matrices
n_params = size(mcmc_params,2)
A2 = zeros(n_samples, n_params)
B2 = zeros(n_samples, n_params)

# sample values for each parameter to create A2 and B2 matrices
for i in 1:n_params
    for j in 1:n_samples
        A2[j,i] = sample_value(param_val=mcmc_params[j,i], bw=bandwidth[i], bounds=bounds[i])
        B2[j,i] = sample_value(param_val=mcmc_params[j,i], bw=bandwidth[i], bounds=bounds[i])
    end
end

#----------------------------------------------- Concatenate design matrices A1 & A2, and B1 & B2 ------------------------------------------------------------------------------------#

# concatenate matrices to get final A and B design matrices
A = hcat(A1', A2) # combine A1 and A2 to get design matrix A
B = hcat(B1', B2) # combine B1 and B2 to get design matrix B

# save parameter names for df
param_names = ["gamma_g", "t_peak", "gamma_d", names(mcmc_params)...]

# save design matrices as df to add parameter names
A_df = DataFrame(A, param_names)
B_df = DataFrame(B, param_names)

#----------------------------------------------- Adjust parameters and save files ----------------------------------------------------------------------------------------------------#

# remove columns for unused parameters (σ²_white_noise_CO₂ and α₀_CO₂)
select!(A_df, Not([:sigma_whitenoise_co2, :alpha0_CO2]))
select!(B_df, Not([:sigma_whitenoise_co2, :alpha0_CO2]))

# add column for land water storage random sample
A_df[!, :lw_random_sample] = rand(Normal(0.0003, 0.00018), n_samples)
B_df[!, :lw_random_sample] = rand(Normal(0.0003, 0.00018), n_samples)

# write A and B design matrices to .csv files (each with size n_samples x n_params)
save(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "sobol_input_A.csv"), A_df)
save(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "sobol_input_B.csv"), B_df)