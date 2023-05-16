using JLACovarianceMatrix
using BetterInputFiles
using Test

@testset "JLACovarianceMatrix.jl" begin
    test_toml = setup_input("Inputs/test.toml", false)
    test_covariance_matrix = main(test_toml)
    
    test_load_toml = setup_input("Inputs/test_load.toml", false)
    test_load_covariance_matrix = main(test_load_toml)

    @test isapprox(test_covariance_matrix, test_load_covariance_matrix)
end
