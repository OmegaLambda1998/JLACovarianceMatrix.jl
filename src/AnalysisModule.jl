# Analysis Module
module AnalysisModule

# Internal Packages 
using ..CovarianceModule

# External Packages 
using ProgressMeter
using KernelDensity
using CairoMakie
CairoMakie.activate!(type="svg")

# Exports
export plot_covariance_matrix

# Constants
const COLOURS::Vector{String} = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"]
const COLOUR_TRUTH::String = "#1b9e77"
const COLOUR_BESTFIT_68::String = "#d95f02"
const COLOUR_BESTFIT_95::String = "#7570b3"
const COLOUR_FIT_68::String = "#e6ab02"
const COLOUR_FIT_95::String = "#a6761d"

function get_uncertainty_sample(covariance_matrix::CovarianceMatrix)
    @info "Generating uncertainty sample"
    distribution = generateDistribution(covariance_matrix)
    rand_draws = rand(distribution, 10000)
    filter = Dict("filter_" * k => rand_draws[i, :] for (i, k) in enumerate(covariance_matrix.keys))
    zp = Dict("zp_" * k => rand_draws[i+length(covariance_matrix.keys), :] for (i, k) in enumerate(covariance_matrix.keys))
    return filter, zp 
end


function plot_covariance_matrix(covariance_matrix::CovarianceMatrix, output::AbstractString)
    filter, zp = get_uncertainty_sample(covariance_matrix)
    uncertainties = Dict("Filter" => filter, "ZP" => zp)
    for (name, uncertainty) in uncertainties
        @info "Plotting covariance matrix"
        fig = Figure(resolution=(3200, 2400), fontsize=28)
        gax = fig[1, 1] = GridLayout()
        n = length(keys(uncertainty))
        ks = sort(collect(keys(uncertainty)))
        p = Progress(n)
        for i in 1:n
            for j in 1:i
                ax = Axis(gax[i, j], xlabel=ks[j], ylabel=ks[i])
                i_chain = uncertainty[ks[i]]
                j_chain = uncertainty[ks[j]]
                if i == j
                    density!(ax, i_chain, color=COLOUR_FIT_68)
                    if i == 1
                        hideydecorations!(ax)
                        if i != n
                            hidexdecorations!(ax, grid=false)
                        end
                    end
                else
                    k = kde(hcat(j_chain, i_chain))
                    x = k.x
                    y = k.y
                    z = k.density ./ maximum(k.density)
                    contour!(ax, x, y, z, levels=[1 - 0.95, 1 - 0.68], colormap=[COLOUR_FIT_68, COLOUR_FIT_95])
                    if j != 1
                        hideydecorations!(ax, grid=false)
                    end
                    if i != n
                        hidexdecorations!(ax, grid=false)
                    end
                end
            end
            next!(p)
        end
        @info "Saving plot"
        output_file = joinpath(output, "$(name)_contour.svg")
        save(output_file, fig)

    end
end

end
