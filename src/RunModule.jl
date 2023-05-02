module RunModule

# External Packages
using Distributions

# Internal Packages
include("InstrumentModule.jl")
using .InstrumentModule

include("CovarianceModule.jl")
using .CovarianceModule

# Exports
export run_JLACovarianceMatrix

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
        analysis = toml["ANALYSIS"]
        if "DRAW" in keys(analysis)
            cov_mat = generateMatrix(covariance_matrix)
            d = MvNormal(zeros(Float64, size(cov_mat, 1)), cov_mat)
            rand_draw = rand(d, analysis["DRAW"])
            filter = Dict(k => rand_draw[i] for (i, k) in enumerate(covariance_matrix.keys))
            zp = Dict(k => rand_draw[i+length(covariance_matrix.keys)] for (i, k) in enumerate(covariance_matrix.keys))
            @show filter
            @show zp
        end
    end
end

end
