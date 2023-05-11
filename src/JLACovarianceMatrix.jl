module JLACovarianceMatrix

# External packages
using TOML
using BetterInputFiles
using ArgParse

# Internal Packages
include("RunModule.jl")
using .RunModule: run_JLACovarianceMatrix

# Exports
export JLACovarianceMatrix

function julia_main()::Cint
    try
        main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function get_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--verbose", "-v"
        help = "Increase level of logging verbosity"
        action = :store_true
        "input"
        help = "Path to .toml file"
        required = true
    end
    return parse_args(s)
end

function JLACovarianceMatrix()
    args = get_args()
    verbose = args["verbose"]
    toml_path = args["input"]
    toml = setup_input(toml_path, verbose)
    return JLACovarianceMatrix(toml)
end

function JLACovarianceMatrix(toml::Dict{String,Any})
    return run_JLACovarianceMatrix(toml)
end

if abspath(PROGRAM_FILE) == @__FILE__
    JLACovarianceMatrix()
end

end
