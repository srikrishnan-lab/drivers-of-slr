using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using GlobalSensitivity
using CSVFiles
using DataFrames

include("gsa_functions.jl")

# define size of ensemble and number of bootstrap samples
#n_samples = 1000 #100000
n_boot = 10 #10000 # one order of magnitude below n_samples # try 100
# nboot needs to be less than the number of COLUMNS in A and B

# read in the formatted design matrices A and B
A = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_input_A.csv")))
B = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_input_B.csv")))
select!(A, Not([:sigma_whitenoise_co2, :alpha0_CO2]))
select!(B, Not([:sigma_whitenoise_co2, :alpha0_CO2]))

# intialize values for years to analyze in sensitivity analysis
gsa_first_yr = 2030 # don't need analyze historical data
gsa_last_yr = 2300

# vector of years that we want to consider
yrs = collect(gsa_first_yr:100:gsa_last_yr)

# -------------------------------------------------------------------------------------------------------------------------------------------- #


# function that takes in an input matrix M where each row is a set of parameters, and then returns a vector of GMSLR for a given year
function brick_run(M::Matrix{Float64}; yr=2100)
    # intialize values
    start_year = 1850
    end_year = 2300

    # function returns a df of global mean sea level rise values
    gmslr = model_ensemble(M, start_year=start_year, end_year=end_year)

    return gmslr # return GMSLR for the specified year (vector with length n_samples)
end

#initialize values
n_params = size(A,2)
n_years = length(yrs)
param_names = names(A)

# first and total order have size (n_params)
first_order = zeros(n_params, n_years)
first_CI = zeros(n_params, n_years)
total_order = zeros(n_params, n_years)
total_CI = zeros(n_params, n_years)
# second order has size (n_params x n_params)
#second_order = zeros(n_params, n_params)
#second_CI = zeros(n_params, n_params)

# conduct global sensitivity analysis
for i in 1:length(yrs)
    
    sobol_output = gsa(M -> brick_run(M; yr=2), Sobol(order=[0,1]), transpose(Matrix(A)[1:100,:]), transpose(Matrix(B)[1:100,:]), batch=true)

    # first and total order
    first_order[:,i] = sobol_output.S1
    first_CI[:,i] = sobol_output.S1_Conf_Int
    total_order[:,i] = sobol_output.ST
    total_CI[:,i] = sobol_output.ST_Conf_Int
    # second order
    #second_order = sobol_output.S2
    #second_CI = sobol_output.S2_Conf_Int
end

# convert Matrix to df and add column names (years)
first_order_df = DataFrame(first_order, Symbol.([yrs...]))
first_CI_df = DataFrame(first_CI, Symbol.([yrs...]))
total_order_df = DataFrame(total_order, Symbol.([yrs...]))
total_CI_df = DataFrame(total_CI, Symbol.([yrs...]))
# second order
#second_order_df = DataFrame(second_order, param_names)
#second_CI_df = DataFrame(second_CI, param_names)

# add a column to each df containing parameter names
insertcols!(first_order_df, 1, :parameter => param_names)
insertcols!(first_CI_df, 1, :parameter => param_names)
insertcols!(total_order_df, 1, :parameter => param_names)
insertcols!(total_CI_df, 1, :parameter => param_names)
# second order
#insertcols!(second_order_df, 1, :parameter => param_names)
#insertcols!(second_CI_df, 1, :parameter => param_names)

# save the results to csv files
save(joinpath(@__DIR__, "..", "results", "sa_results", "test", "first_order.csv"), first_order_df)
save(joinpath(@__DIR__, "..", "results", "sa_results", "test", "first_CI.csv"), first_CI_df)
save(joinpath(@__DIR__, "..", "results", "sa_results", "test", "total_order.csv"), total_order_df)
save(joinpath(@__DIR__, "..", "results", "sa_results", "test", "total_CI.csv"), total_CI_df)
# second order
#save(joinpath(@__DIR__, "..", "results", "sa_results", "test", "second_order", "sens_2150.csv"), second_order_df)
#save(joinpath(@__DIR__, "..", "results", "sa_results", "test", "second_order", "CI_2150.csv"), second_CI_df)
