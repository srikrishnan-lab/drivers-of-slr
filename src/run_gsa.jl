using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using GlobalSensitivity

include("gsa_functions.jl")

# define size of ensemble and number of bootstrap samples
n_samples = 100 #100000
n_boot = 100 #10000 # one order of magnitude below n_samples

# intialize values for years
current_year = 2030 # first year for sensitivity analysis (don't need analyze historical data)
end_year = 2300

# read in the design matrices A and B
A = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_input_A.csv")))
B = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_input_B.csv")))

# -------------------------------------------------------------------------------------------------------------------------------------------- #
#attempt for batch=false

# function that takes in an input vector V that is a set of parameters, and then returns a vector of GMSLR for a given year
function f_false(V)
    # intialize values
    start_year = 1850
    end_year = 2300

    # function returns a df of global mean sea level rise values from 1850-2300
    gmslr = model_ensemble(num_samples=1, calibrated_params=V)

    # subset df to get values for desired year
    col = findfirst(2100 .∈ start_year:end_year) # gets index for column of specified year
    yr_gmslr = gmslr[:,col]
       
    return yr_gmslr # return GMSLR for the specified year (vector with length n_samples)
end

# batch=false test (this works)
my_gmslr_2 = model_ensemble(num_samples=1, calibrated_params=Vector(A[1,:]))
my_yr_2 = f_false(Vector(A[1,:]))

# conduct global sensitivity analysis
gsa(f_false, Sobol(order=[0,1], nboot=n_boot), Matrix(A), Matrix(B), batch=false)

# -------------------------------------------------------------------------------------------------------------------------------------------- #
# attempt for batch=true
# function that takes in an input matrix M where each row is a set of parameters, and then returns a vector of GMSLR for a given year
function f(M; n_samples, year=2100)
    # intialize values
    start_year = 1850
    end_year = 2300

    # function returns a df of global mean sea level rise values from 1850-2300
    gmslr = model_ensemble(num_samples=n_samples, calibrated_params=M)

    # subset df to get values for desired year
    col = findfirst(year .∈ start_year:end_year) # gets index for column of specified year
    yr_gmslr = gmslr[:,col]
       
    return yr_gmslr # return GMSLR for the specified year (vector with length n_samples)
end

# vector of years that we want to consider
yrs = collect(current_year:10:end_year)

# batch=true test (this works)
my_gmslr = model_ensemble(num_samples=n_samples, calibrated_params=Matrix(A))
my_yr = f(Matrix(A), n_samples=n_samples)

# conduct global sensitivity analysis
gsa(M -> f(M, n_samples=n_samples, year=2100), Sobol(order=[0,1], nboot=n_boot), Matrix(A), Matrix(B), batch=true)

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# batch=true test (this works)
my_gmslr = model_ensemble(num_samples=n_samples, calibrated_params=Matrix(A))
my_yr = f(Matrix(A), n_samples=n_samples)
col = findfirst(2100 .∈ 1850:2300)
my_gmslr[:, col]

function new_fun(M; year=2100)
    for row in eachrow(M)
        model_ensemble(num_samples=n_samples, calibrated_params=M)
        spit out sea level
    end
    return gmslr_out (vector of n_samples length)
end
# vectors of years every 5-10 years from 1850-2300
sobol[:,i] = gsa(M -> f(M, yr=yr[i]), Sobol(order=[0,1], nboot=n_boot), A, B, batch=true) # this will output n_years x n_parameters

gsa = gsa(f, Sobol(order=[0,1], nboot=n_boot), A, B, batch=false) # this will output n_years x n_parameters

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# below this is old stuff, not relevant
#=
# function f(M) that takes in a matrix M, which is the same size as A and B
# A and B and M are all size (n_samples x n_years)
function old_f(M; A=A, B=B, year=2100)
    output = zeros(dim(M)[1]) # output is out_M
    C = vcat(A,B)
    C_out = vcat(A_out, B_out)
    for i=1:dim(M)[1] # size(M,1) each row in M
        j=find(M[i] in C) # look up index in C, only need to find first (not all) -> find index where row in M matches a row in C
        out[i] = C_out[j]
    end
    return out
end

#=
for i in 1:size(M,1) # loop through rows of M
    col = findfirst(yr -> yr == year, start_year:end_year) # find index for column containing the year specified
    j = findfirst(val -> val == M[i,col], C_output[:,col]) # find index j where the row in M matches a row in C value in C_output matches M
    findfirst((in)(2010:end_year), start_year:end_year)
end
=#

# read in relevant output (global mean sea level rise) for both A and B
A_gmslr = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_output_A", "gmslr.csv")))
B_gmslr = DataFrame(load(joinpath(@__DIR__, "..", "results", "sa_results", "sobol_output_B", "gmslr.csv")))

#start_year = 1850
#end_year = 2300

function old_f(M_in; A_in=A_base, B_in=B_base, A_out=A_gmslr, B_out=B_gmslr, yr=2100) # takes in M_in, which is a matrix of parameters (same size as A_in and B_in)

    # concatenate A and B parameter inputs to get matrix C_in (n_samples x n_params)
    C_in = hcat(A_in', B_in')
    # concatenate A and B GMSLR output to get matrix C_out (n_samples x n_years)
    C_out = vcat(A_out, B_out)

    # initialize matrix to store output (n_samples x n_years)
#    M_out = zeros(size(A_out,1), length(start_year:end_year))

    # function that returns the index of the first row in the df (C_in) that matches the inputted row (in M_in)
    match_row(row, df) = findfirst(i -> all(j -> row[j] == df[i,j], 1:size(df,2)), 1:size(df,1))

    # loop through rows of M_in (each row is a set of samples)
#    n_rows = size(M_in, 1)
 #   for i in 1:n_rows
#        j = match_row(M_in[i,:], C_in) # j is the desired row index
        j = match_row(M_in, C_in') # j is the desired row index

        col = findfirst(yr .∈ start_year:end_year) # find index for column containing the year specified
        return C_out[j,col]
#    end
    # return the output (n_samples x n_years)
    return M_out
end
=#