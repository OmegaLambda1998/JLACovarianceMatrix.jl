# Covariance Module
module CovarianceModule

# Internal Packages 
using ..InstrumentModule

# External Packages 
using JLD2
using LinearAlgebra
using Combinatorics
using Distributions

# Exports
export CovarianceMatrix
export saveCovarianceMatrix
export loadCovarianceMatrix
export generateMatrix
export generateDistribution
export draw_covariance_matrix
export get_uncertainty_sample

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

struct CovarianceMatrix
    keys::Vector{String}
    ZPUncertainty::Matrix{Float64} # Default units of Mag^2 
    FilterUncertainty::Vector{Float64}
end

function Base.isapprox(cov1::CovarianceMatrix, cov2::CovarianceMatrix; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps(), nans::Bool=false)
    return isapprox(cov1.ZPUncertainty, cov2.ZPUncertainty; atol=atol, rtol=rtol, nans=nans) && isapprox(cov1.FilterUncertainty, cov2.FilterUncertainty; atol=atol, rtol=rtol, nans=nans)
end


# Constants
const SLOPE = 0.005
const WAVE_START = 300.0 # nm
const WAVE_END = 1000.0 # nm
const CENTRAL = 555.6 # nm
const VAL = SLOPE / (WAVE_END - WAVE_START)
const SUPERCAL_UNCERTAINTY = 1.0 / 3.0

function CovarianceMatrix(instruments::Dict{String,Instrument})
    filters = Dict{String,Filter}()
    for instrument in values(instruments)
        merge!(filters, instrument.filters)
    end
    ndim = length(filters)
    keys = Vector{String}(undef, ndim)
    ZPUncertainty = zeros(Float64, ndim, ndim)
    FilterUncertainty = Vector{Float64}(undef, ndim)
    for (i, (filter_name, filter_i)) in enumerate(filters)
        keys[i] = filter_name
        SuperCalCorrection = (((1 + SUPERCAL_UNCERTAINTY) * filter_i.SuperCalCorrection)^2)
        # Filter Uncertainty already squared
        FilterUncertainty[i] = filter_i.FilterUncertainty + SuperCalCorrection
        for (j, (filter_name_j, filter_j)) in enumerate(filters)
            if i == j
                # Already squared
                ZPUncertainty[i, j] = filter_i.ZPUncertainty
            end
            diff_i = filter_i.CentralWavelength - CENTRAL
            diff_j = filter_j.CentralWavelength - CENTRAL
            ZPUncertainty[i, j] += VAL * VAL * diff_i * diff_j
        end
    end
    return CovarianceMatrix(keys, ZPUncertainty, FilterUncertainty)
end

function saveCovarianceMatrix(covariance::CovarianceMatrix, output::AbstractString)
    save(output, Dict(string(key) => getfield(covariance, key) for key in fieldnames(CovarianceMatrix)))
end

function loadCovarianceMatrix(input::AbstractString)
    d = load(input)
    val = [d[string(key)] for key in fieldnames(CovarianceMatrix)]
    return CovarianceMatrix(val...)
end

function generateMatrix(covariance::CovarianceMatrix)
    nfilt = length(covariance.FilterUncertainty)
    FilterUncertaintyMatrix = diagm(covariance.FilterUncertainty)
    ZeroMatrix = zeros(Float64, nfilt, nfilt)
    cov_matrix = [FilterUncertaintyMatrix ZeroMatrix; ZeroMatrix covariance.ZPUncertainty]
    if !(isequal(cov_matrix, transpose(cov_matrix)))
        if isapprox(cov_matrix, transpose(cov_matrix))
            @info "Matrix is not symmetric, but is approximately symmetric, fixing."
            cov_matrix = 0.5 * (cov_matrix + transpose(cov_matrix))
        else
            @info "Matrix is not symmetric"
            non_symmetric_indices = findall(cov_matrix .!= transpose(cov_matrix))
            for (i, j) in Tuple.(non_symmetric_indices)
                @debug "i: $i, j: j"
                @debug "Covariance[i,j]: $(cov_matrix[i, j])"
                @debug "Covariance[j,i]: $(cov_matrix[j, i])"
            end
        end
    end
    eigen_values = eigvals(cov_matrix)
    if any(eigen_values .< 0)
        @info "Matrix is not positive definite, fixing."
        offset_matrix = I * (eps(Float64) - minimum(eigen_values))
        @debug "$offset_matrix"
        cov_matrix += offset_matrix
    end
    return cov_matrix
end

function generateDistribution(covariance::CovarianceMatrix)
    cov_mat = generateMatrix(covariance)
    distribution = MvNormal(zeros(Float64, size(cov_mat, 1)), cov_mat)
    return distribution
end

function get_uncertainty_sample(covariance_matrix::CovarianceMatrix, num=10000)
    @info "Generating uncertainty sample"
    distribution = generateDistribution(covariance_matrix)
    rand_draws = rand(distribution, num)
    filter = Dict(k => rand_draws[i, :] for (i, k) in enumerate(covariance_matrix.keys))
    zp = Dict(k => rand_draws[i+length(covariance_matrix.keys), :] for (i, k) in enumerate(covariance_matrix.keys))
    return filter, zp
end

function draw_covariance_matrix(covariance_matrix::CovarianceMatrix, config::Dict{String,Any}; output::Union{AbstractString,Nothing}=nothing)
    num = get(config, "NUM", 100)
    SALTJacobian = get(config, "SALTJACOBIAN", nothing)
    d_filter, d_zp = get_uncertainty_sample(covariance_matrix, num)
    if !isnothing(SALTJacobian)
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
            band = spl[end]
            instrument = MAPPING[base_key]
            mag = get(samemag, instrument, instrument)
            filter = get(samefilter, instrument, instrument)
            d_filter["$(filter)_$(band)"] = d_filter[key]
            delete!(d_filter, key)
            d_zp["$(mag)_$(band)"] = d_zp[key]
            delete!(d_zp, key)
        end
    end
    if !isnothing(output)
        d_filter_file = joinpath(output, "FilterUncertainty.jld2")
        @info "Saving filter uncertainties to $d_filter_file"
        d_zp_file = joinpath(output, "ZPUncertainty.jld2")
        @info "Saving zp uncertainties to $d_zp_file"
        save(d_filter_file, d_filter)
        save(d_zp_file, d_zp)
    end
    return d_filter, d_zp
end

end
