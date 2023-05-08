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


struct CovarianceMatrix
    keys::Vector{String}
    ZPUncertainty::Matrix{Float64} # Default units of Mag^2 
    FilterUncertainty::Vector{Float64}
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

end
