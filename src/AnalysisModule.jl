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
export draw_covariance_matrix

# Constants
const COLOURS::Vector{String} = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"]
const COLOUR_TRUTH::String = "#1b9e77"
const COLOUR_BESTFIT_68::String = "#d95f02"
const COLOUR_BESTFIT_95::String = "#7570b3"
const COLOUR_FIT_68::String = "#e6ab02"
const COLOUR_FIT_95::String = "#a6761d"

# Convert from CovarianceMatrix name to SALT2 name (which then gets converted to SNANA name later on)
const MAPPING = Dict{String,String}(
    "CFA3K" => "CfA3_KEPLERCAM",
    "CFA3S" => "CfA3_STANDARD",
    "CFA4_1" => "CfA1",
    "CFA4_2" => "CfA2",
    "CSP" => "CSP",
    "DES" => "DES",
    "SDSS" => "SDSS",
    "SNLS" => "SNLS",
    "STANDARD" => "STANDARD"
)

function get_uncertainty_sample(covariance_matrix::CovarianceMatrix, num=10000)
    @info "Generating uncertainty sample"
    distribution = generateDistribution(covariance_matrix)
    rand_draws = rand(distribution, num)
    filter = Dict(k => rand_draws[i, :] for (i, k) in enumerate(covariance_matrix.keys))
    zp = Dict(k => rand_draws[i+length(covariance_matrix.keys), :] for (i, k) in enumerate(covariance_matrix.keys))
    return filter, zp
end


function plot_covariance_matrix(covariance_matrix::CovarianceMatrix, config::Dict{String,Any}, output::AbstractString)
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

function draw_covariance_matrix(covariance_matrix::CovarianceMatrix, config::Dict{String,Any}, output::AbstractString)
    num = get(config, "NUM", 100)
    SALTJacobian = get(config, "SALTJACOBIAN", nothing)
    filter, zp = get_uncertainty_sample(covariance_matrix, num)
    if isnothing(SALTJacobian)
        @info "Filter uncertainties:\n$filter"
        @info "ZP uncertainties:\n$zp"
    else
        samemag_list = get(config, "SURVEY_LIST_SAMEMAGSYS", Vector{String}())
        if length(samemag_list) > 0
            samemag = Dict{String,String}(mag => samemag_list[1] for mag in samemag_list[2:end])
        else
            samemag = Dict{String,String}()
        end
        samefilter_list = get(config, "SURVEY_LIST_SAMEFILTER", Vector{String}())
        if length(samefilter_list) > 0
            samefilter = Dict{String,String}(filter => samefilter_list[1] for filter in samefilter_list[2:end])
        else
            samefilter = Dict{String,String}()
        end
        for key in sort(covariance_matrix.keys)
            spl = split(key, "_")
            base_key = join(spl[1:end-1], "_")
            @show base_key
            band = spl[end]
            instrument = MAPPING[base_key]
            mag = get(samemag, instrument, instrument)
            filter = get(samefilter, instrument, instrument)
            @show mag, filter, band
        end
    end
end

end
