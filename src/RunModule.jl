module RunModule

# External Packages

# Internal Packages
include(joinpath(@__DIR__, "InstrumentModule.jl"))
using .InstrumentModule

include(joinpath(@__DIR__, "CovarianceModule.jl"))
using .CovarianceModule



# Exports
export run_JLACovarianceMatrix
export CovarianceMatrix
export loadCovarianceMatrix, saveCovarianceMatrix
export draw_covariance_matrix

function run_JLACovarianceMatrix(toml::Dict{String,Any})
    config = get(toml, "COVARIANCEMATRIX", Dict{String,Any}())
    name = get(config, "NAME", "DES")
    if "INSTRUMENT" in keys(config)
        instruments = Dict{String,Instrument}()
        @info "Loading in instruments"
        for (instrument_name, instrument_dict) in config["INSTRUMENT"]
            @info "Loading in $instrument_name"
            instruments[instrument_name] = Instrument(instrument_name, instrument_dict)
        end
        @info "Finished loading in instruments"
        @info "Creating Covariance Matrix"
        covariance_matrix = CovarianceMatrix(instruments)
        @info "Finished creating Covariance Matrix"
    elseif "INPUT" in keys(config)
        input = config["INPUT"]
        if !isabspath(input)
            input = joinpath(toml["GLOBAL"]["BASE_PATH"], config["INPUT"])
        end
        input = abspath(input)
        @info "Loading Covariance Matrix from $input"
        covariance_matrix = loadCovarianceMatrix(input)
        @info "Finished loading Covariance Matrix"
    else
        error("Please either specify a premade covariance matrix via INPUT, or create a new covariance matrix via INSTRUMENT")
    end
    saveCovarianceMatrix(covariance_matrix, joinpath(toml["GLOBAL"]["OUTPUT_PATH"], "$(name).jld2"))
    if "ANALYSIS" in keys(toml)
        for analysis in toml["ANALYSIS"]
            output = get(analysis, "OUTPUT", "Output")
            if !isabspath(output)
                output = joinpath(toml["GLOBAL"]["OUTPUT_PATH"], output)
            end
            output = abspath(output)
            if !isdir(output)
                mkdir(output)
            end
            if "PLOT" in keys(analysis)
                # Only include PlotModule if needed to save on precompilation time if we don't need to use Makie
                include(joinpath(@__DIR__, "PlotModule.jl"))
                @eval using .PlotModule

                @info "Plotting Covariance Matrix"
                Base.invokelatest(plot_covariance_matrix, covariance_matrix, analysis["PLOT"], output)
            end
            if "DRAW" in keys(analysis)
                @info "Drawing from Covariance Matrix"
                d_filter, d_zp = draw_covariance_matrix(covariance_matrix, analysis["DRAW"]; output=output)
            end
        end
    end
    return covariance_matrix
end

end
