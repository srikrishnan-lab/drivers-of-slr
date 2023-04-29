using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using GlobalSensitivity

include("gsa_functions.jl")

# define size of ensemble and number of bootstrap samples
n_samples = 100 #100000
n_boot = 100 #10000 # one order of magnitude below n_samples

# read in the design matrices A and B
A = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_input_A.csv")))
B = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_input_B.csv")))

# intialize values for years to analyze in sensitivity analysis
gsa_first_yr = 2030 # don't need analyze historical data
gsa_last_yr = 2300

# vector of years that we want to consider
yrs = collect(2100:100:2300)

# specify directory to store results
output_dir = "test"

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

#initialize values
n_params = size(A,2)
n_years = length(yrs)
param_names = names(A)

# first and total order have size (n_params) for each year
first_order = zeros(n_params, n_years)
first_CI = zeros(n_params, n_years)
total_order = zeros(n_params, n_years)
total_CI = zeros(n_params, n_years)
# second order has size (n_params x n_params) for each year

# loop through years we want to analyze
for i in 1:length(yrs)

    # conduct global sensitivity analysis
    sobol_output = gsa(M -> brick_run(M; yr=yrs[i]), Sobol(order=[0,1,2], nboot=2), transpose(Matrix(A)[1:5,:]), transpose(Matrix(B)[1:5,:]), batch=true)

    # store first and total order sensitivities and confidence intervals
    first_order[:,i] = sobol_output.S1
    first_CI[:,i] = sobol_output.S1_Conf_Int
    total_order[:,i] = sobol_output.ST
    total_CI[:,i] = sobol_output.ST_Conf_Int

    # store second order sensitivities and confidence intervals
    second_order = sobol_output.S2
    second_CI = sobol_output.S2_Conf_Int
    # convert Matrix to df and add column names (parameter names)
    second_order_df = DataFrame(second_order[:,:,1], param_names)
    second_CI_df = DataFrame(second_CI[:,:,1], param_names)
    # add a column to each df containing parameter names
    insertcols!(second_order_df, 1, :parameter => param_names)
    insertcols!(second_CI_df, 1, :parameter => param_names)
    # save the results to csv files
    current_yr = yrs[i]
    save(joinpath(@__DIR__, "..", "results", "sa_results", "$output_dir", "second_order", "sens_$current_yr.csv"), second_order_df)
    save(joinpath(@__DIR__, "..", "results", "sa_results", "$output_dir", "second_order", "CI_$current_yr.csv"), second_CI_df)

end

# convert Matrix to df and add column names (years)
first_order_df = DataFrame(first_order, Symbol.([yrs...]))
first_CI_df = DataFrame(first_CI, Symbol.([yrs...]))
total_order_df = DataFrame(total_order, Symbol.([yrs...]))
total_CI_df = DataFrame(total_CI, Symbol.([yrs...]))

# add a column to each df containing parameter names
insertcols!(first_order_df, 1, :parameter => param_names)
insertcols!(first_CI_df, 1, :parameter => param_names)
insertcols!(total_order_df, 1, :parameter => param_names)
insertcols!(total_CI_df, 1, :parameter => param_names)

# save the results to csv files
save(joinpath(@__DIR__, "..", "results", "sa_results", "$output_dir", "first_order.csv"), first_order_df)
save(joinpath(@__DIR__, "..", "results", "sa_results", "$output_dir", "first_CI.csv"), first_CI_df)
save(joinpath(@__DIR__, "..", "results", "sa_results", "$output_dir", "total_order.csv"), total_order_df)
save(joinpath(@__DIR__, "..", "results", "sa_results", "$output_dir", "total_CI.csv"), total_CI_df)