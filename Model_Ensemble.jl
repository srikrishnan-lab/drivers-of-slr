# activate the environment
using Pkg
Pkg.activate(joinpath(@__DIR__, "sea-level-rise-env/")) # this does the same thing: Pkg.activate("/Users/ced227/Documents/Repositories/drivers-of-sea-level-rise/sea-level-rise-env/")
Pkg.instantiate()

# import necessary packages
using Mimi
using MimiDICE2013
using MimiSNEASY
using MimiBRICK
using DataFrames
using Impute

print(T)

# run DICE
DICE = MimiDICE2013.get_model()
run(DICE)
#explore(DICE)

# isolate output from DICE
co2_emissions = getdataframe(DICE, :emissions=>:E) # E is Total CO2 emissions (GtCO2 per year) (5 year timesteps)

# create new df to fill in missing values for 1-year timesteps
interpolated_df = DataFrame(Year = 2010:2305, Emissions = missing)
# convert Emissions column type to both Missing and Float so that Float values can replace missing values
interpolated_df[!, :Emissions] = convert(Vector{Union{Float64,Missing}}, interpolated_df[:, :Emissions])

# loop through DICE df output
for i in 1:nrow(co2_emissions)
    emission_value = co2_emissions[i, :E]
    loc = findfirst(==(co2_emissions[i, :time]), interpolated_df[:, :Year])
    interpolated_df[loc,:Emissions] = emission_value
end

# interpolate to fill in missing values in interpolated_df
final_co2_df = Impute.interp(interpolated_df)
final_emissions = final_co2_df[: , 2]

# create SNEASY-BRICK
SNEASY_BRICK = MimiBRICK.create_sneasy_brick(start_year=2010, end_year=2305)

# update parameters to feed DICE emissions to SNEASY-BRICK
update_param!(SNEASY_BRICK, :ccm, :CO2_emissions, final_emissions)
run(SNEASY_BRICK)
#explore(SNEASY_BRICK)

# retrieve output from SNEASY-BRICK
sea_level_rise = getdataframe(SNEASY_BRICK, :global_sea_level=>:sea_level_rise) # total sea level rise from all components (m) (includes landwater storage for projection periods)

# final output dataframe
output_df = sea_level_rise
insertcols!(output_df, 2, :co2_emissions => final_emissions)
CSV.write("/home/fs01/ced227/drivers-of-sea-level-rise/output.csv", output_df,  header = false)
