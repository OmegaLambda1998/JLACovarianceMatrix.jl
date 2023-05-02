module InstrumentModule

# External Packages

# Internal Packages

# Exports
export Instrument
export Filter

struct Filter
    name::String
    ZPUncertainty::Float64 # Default unit of mag^2
    SuperCalCorrection::Float64 # Default unit of nm
    FilterUncertainty::Float64 # Default unit of nm
    CentralWavelength::Float64 # Default unit of nm
end


struct Instrument
    filters::Dict{String,Filter}
end

function Instrument(instrument_name::String, instrument_dict::Dict{String,Any})
    filters = Dict{String,Filter}()
    filter_names = instrument_dict["FILTERS"]
    if "ZPUNCERTAINTY" in keys(instrument_dict)
        ZPUncertainties = instrument_dict["ZPUNCERTAINTY"] .^ 2
    else
        astrometry = instrument_dict["ASTROMETRY"]
        nonlinearity = instrument_dict["NONLINEARITY"]
        photometric_zero_pointing = instrument_dict["PHOTOMETRICZEROPOINTING"]
        photometric_bias = instrument_dict["PHOTOMETRICBIAS"]
        uniformity = instrument_dict["UNIFORMITY"]
        base_uncertainty = (astrometry * astrometry) + (nonlinearity * nonlinearity) + (photometric_zero_pointing * photometric_zero_pointing) + (photometric_bias * photometric_bias) + (uniformity * uniformity)
        ABUncertainties = instrument_dict["ABUNCERTAINTY"] .^ 2
        ZPUncertainties = ((base_uncertainty * base_uncertainty) .+ ABUncertainties)
    end
    SuperCalCorrections = instrument_dict["SUPERCALCORRECTION"]
    FilterUncetainties = instrument_dict["FILTERUNCERTAINTY"] .^ 2
    CentralWavelengths = instrument_dict["CENTRALWAVELENGTH"]
    for (i, name) in enumerate(filter_names)
        f_name = instrument_name * "_" * name
        filters[f_name] = Filter(f_name, ZPUncertainties[i], SuperCalCorrections[i], FilterUncetainties[i], CentralWavelengths[i])
    end
    return Instrument(filters)
end

end
