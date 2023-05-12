# this script plots the results from the Global Sensitivity Analysis

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Plots, Measures
include("../src/functions.jl")

output_dir = "10_000_samples"

# read in results
first_order = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "first_order.csv")))
first_CI    = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "first_CI.csv")))
total_order = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "total_order.csv")))
total_CI    = DataFrame(load(joinpath(@__DIR__, "..", "results", "gsa_results", "$output_dir", "total_CI.csv")))

# treat negative sensitivities as zero
first_order[:,2:end] .= ifelse.(first_order[:,2:end] .< 0, 0, first_order[:,2:end])
total_order[:,2:end] .= ifelse.(total_order[:,2:end] .< 0, 0, total_order[:,2:end])

# -------------------------------------- Plot first order sensitivities --------------------------------------- #
all_plts = []
years = ["2100", "2200", "2300"]

for year in years
    plt = scatter(first_order.parameter, first_order[:,"$year"], yerror = first_CI[:,"$year"], 
            title = "$year", xlabel="Parameter", ylabel="Sensitivity", xticks=:all, legend=:false,
            xrotation=45, markersize=3, tickfontsize=6, bottom_margin=7mm, left_margin=5mm)
    push!(all_plts, plt)
end

p = plot(all_plts..., plot_title="First Order Sensitivities", layout=(3,1), size=(900,1200))
#savefig(p, "/Users/ced227/Desktop/plots/first_order.pdf")

# -------------------------------------- Plot total order sensitivities --------------------------------------- #
all_plts = []
years = ["2100", "2200", "2300"]

for year in years
    plt = scatter(total_order.parameter, total_order[:,"$year"], yerror = total_CI[:,"$year"], 
            title = "$year", xlabel="Parameter", ylabel="Sensitivity", xticks=:all, legend=:false,
            xrotation=45, markersize=3, tickfontsize=6, bottom_margin=7mm, left_margin=5mm)
    push!(all_plts, plt)
end

plot(all_plts..., plot_title="Total Order Sensitivities", layout=(3,1), size=(900,1200))
#savefig(p, "/Users/ced227/Desktop/plots/total_order.pdf")