using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using GlobalSensitivity
using CSVFiles
using DataFrames

include("gsa_functions.jl")

# define size of ensemble and number of bootstrap samples
#n_samples = 1000 #100000
n_boot = 10 #10000 # one order of magnitude below n_samples 
# nboot needs to be less than the number of COLUMNS in A and B

# read in the formatted design matrices A and B
A = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_input_A.csv")))
B = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_input_B.csv")))

# intialize values for years to analyze in sensitivity analysis
gsa_first_yr = 2030 # don't need analyze historical data
gsa_last_yr = 2300

# vector of years that we want to consider
yrs = collect(gsa_first_yr:10:gsa_last_yr) # yr = yrs[i]

# -------------------------------------------------------------------------------------------------------------------------------------------- #

#function that takes in an input vector v that is a set of parameters, and then returns a vector of GMSLR for a given year
function f(v; yr=2100)
    # intialize values
    start_year = 1850
    end_year = 2300

    # function returns a df of global mean sea level rise values from 1850-2300
    gmslr = model_ensemble(calibrated_params=v, start_year=start_year, end_year=end_year)

    # subset df to get values for desired year
    col = findfirst(yr .∈ start_year:end_year) # gets index for column of specified year
    yr_gmslr = gmslr[:,col]
       
    return yr_gmslr # return GMSLR for the specified year (Float64 scalar)
end

# conduct global sensitivity analysis
sobol_output = gsa(v -> f(v; yr=2100), Sobol(order=[0,1], nboot=2), transpose(Matrix(A)[1:5, :]), transpose(Matrix(B)[1:5, :]))
println(sobol_output)
# -------------------------------------------------------------------------------------------------------------------------------------------- #
# save output
#output = Dict("2100" => Dict("first_order" => sobol_output.S1, "first_CI" =>sobol_output.S1_Conf_Int, "total_order" => sobol_output.ST, "total_CI" => sobol_output.ST_Conf_Int))

#initialize values
n_params = size(A,2)
n_years = 1
param_names = names(A)

# first and total order have size (n_params)
first_order = zeros(n_params, n_years)
first_CI = zeros(n_params, n_years)
total_order = zeros(n_params, n_years)
total_CI = zeros(n_params, n_years)
# second order has size (n_params x n_params)
#second_order = zeros(n_params, n_params)
#second_CI = zeros(n_params, n_params)

# first and total order
first_order[:,1] = sobol_output.S1
first_CI[:,1] = sobol_output.S1_Conf_Int
total_order[:,1] = sobol_output.ST
total_CI[:,1] = sobol_output.ST_Conf_Int
# second order
#second_order = sobol_output.S2
#second_CI = sobol_output.S2_Conf_Int

# convert Matrix to df and add column names (years)
first_order_df = DataFrame(first_order, ["2100"])
first_CI_df = DataFrame(first_CI, ["2100"])
total_order_df = DataFrame(total_order, ["2100"])
total_CI_df = DataFrame(total_CI, ["2100"])
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


#= EVERYTHING BELOW HERE IS OLD AND FOR REFERENCE/IN CASE ONLY
# -------------------------------------------------------------------------------------------------------------------------------------------- #
# !!!!!attempt for batch=TRUE (NOT using matchrow() approach)!!!!!

# function that takes in an input matrix M where each row is a set of parameters, and then returns a vector of GMSLR for a given year
function f(M; yr=2100)
    # intialize values
    start_year = 1850
    end_year = 2300

    # function returns a df of global mean sea level rise values
    gmslr = model_ensemble(calibrated_params=M, start_year=start_year, end_year=end_year)

    # subset df to get values for desired year
    col = findfirst(yr .∈ start_year:end_year) # gets index for column of specified year
    yr_gmslr = gmslr[:,col]
       
    return yr_gmslr' # return GMSLR for the specified year (vector with length n_samples)
end

# conduct global sensitivity analysis
sobol_output = gsa(M -> f(M; yr=2100), Sobol(order=[0,1]), transpose(Matrix(A)[1:10,:]), transpose(Matrix(B)[1:10,:]), batch=true)

# -------------------------------------------------------------------------------------------------------------------------------------------- #
# !!!!!attempt for batch=FALSE (NOT using matchrow() approach)!!!!!

#function that takes in an input vector v that is a set of parameters, and then returns a vector of GMSLR for a given year
function f(v; yr=2100)
    # intialize values
    start_year = 1850
    end_year = 2300

    # function returns a df of global mean sea level rise values from 1850-2300
    gmslr = model_ensemble(calibrated_params=v, start_year=start_year, end_year=end_year)

    # subset df to get values for desired year
    col = findfirst(yr .∈ start_year:end_year) # gets index for column of specified year
    yr_gmslr = gmslr[:,col]
       
    return yr_gmslr # return GMSLR for the specified year (Float64 scalar)
end

# conduct global sensitivity analysis
sobol_output = gsa(v -> f(v; yr=2100), Sobol(order=[0,1], nboot=10), transpose(Matrix(A)[1:10, :]), transpose(Matrix(B)[1:10, :]))

# batch=false test (this works)
my_gmslr_2 = model_ensemble(calibrated_params=Vector(A[1,:]))
my_yr_2 = f(Vector(A[1,:]))

# -------------------------------------------------------------------------------------------------------------------------------------------- #
# !!!!!!! matchrow() approach function WITH transposing!!!!!!!!!

# run function to return a collection of csv files containing corresponding outputs for inputs of design matrices A and B
#model_ensemble(num_samples=n_samples, calibrated_params="sobol_input_A", output_dir="sobol_output_A")
#model_ensemble(num_samples=n_samples, calibrated_params="sobol_input_B", output_dir="sobol_output_B")

# read in relevant output (global mean sea level rise) for both A and B
A_gmslr = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_output_A", "gmslr.csv")))
B_gmslr = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_output_B", "gmslr.csv")))

A_gmslr = transpose(Matrix(A_gmslr))
B_gmslr = transpose(Matrix(B_gmslr))

function transpose_f(M_in; A_in=A_params, B_in=B_params, A_out=A_gmslr, B_out=B_gmslr, year=2100) # takes in M_in, which is a matrix of parameters (same size as A_in and B_in)
    # initialize model years
    start_year = 1850
    end_year = 2300

    # concatenate A and B parameter inputs to get matrix C_in (n_samples x n_params)
    C_in = hcat(A_in, B_in)
    # concatenate A and B GMSLR output to get matrix C_out (n_samples x n_years)
    C_out = hcat(A_out, B_out)

    # initialize vector to store output (length of n_samples)
    output = zeros(size(A_out,2))

    # function that returns the index of the first row in the df (C_in) that matches the inputted row (in M_in)
    match_row(row, df) = findfirst(i -> all(j -> row[j] == df[i,j], 1:size(df,2)), 1:size(df,1))
    match_col(col, df) = findfirst(i -> all(j -> col[j] == df[j,i], 1:size(df,1)), 1:size(df,2))

    # loop through cols of M_in (each column is a set of samples)
    n_cols = size(M_in,2)
    for i in 1:n_cols
        j = match_col(M_in[:,i], C_in) # j is the desired column index
        #j = match_row(M_in, C_in') # j is the desired row index
        row = findfirst(year .∈ start_year:end_year) # find index for column containing the year specified
        output[i] = C_out[row,j]
    end

    # return the output for the year specified (vector with length of n_samples)
    return output'
end

gsa(transpose_f, Sobol(), A_params, B_params, batch=true)

# -------------------------------------------------------------------------------------------------------------------------------------------- #
# !!!!!!! matchrow() approach function WITHOUT transposing!!!!!!!!!

this is right, except for the issue of tranposing. so it's right if gsa didn't need tranpose
function f_no_transpose(M_in; A_in=A_params, B_in=B_params, A_out=A_gmslr, B_out=B_gmslr, year=2100) # takes in M_in, which is a matrix of parameters (same size as A_in and B_in)
    # initialize model years
    start_year = 1850
    end_year = 2300

    # concatenate A and B parameter inputs to get matrix C_in (n_samples x n_params)
    C_in = vcat(A_in, B_in)
    # concatenate A and B GMSLR output to get matrix C_out (n_samples x n_years)
    C_out = vcat(A_out, B_out)

    # initialize vector to store output (length of n_samples)
    output = zeros(size(A_out,1))

    # function that returns the index of the first row in the df (C_in) that matches the inputted row (in M_in)
    match_row(row, df) = findfirst(i -> all(j -> row[j] == df[i,j], 1:size(df,2)), 1:size(df,1))

    # loop through rows of M_in (each row is a set of samples)
    n_rows = size(M_in, 1)
    for i in 1:n_rows
        j = match_row(M_in[i,:], C_in) # j is the desired row index
        #j = match_row(M_in, C_in') # j is the desired row index
        col = findfirst(year .∈ start_year:end_year) # find index for column containing the year specified
        output[i] = C_out[j,col]
    end

    # return the output for the year specified (vector with length of n_samples)
    return output
end

# -------------------------------------------------------------------------------------------------------------------------------------------- #
=#