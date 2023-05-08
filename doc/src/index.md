# JLACovarianceMatrix Documentation

[JLACovarianceMatrix.jl](https://github.com/dessn/JLACovarianceMatrix.jl.git)

Create JLA-like covariance matrices from a `.toml` file. Just provide ZP and Filter uncertainties, alongside (optional) supercal corrections. 

## Install
```julia
using Pkg
Pkg.add("JLACovarianceMatrix")
```

## Usage
Create a new covariance matrix by providing the various uncertainties that make up the total uncertainty budget. Note that the toml keys are case-insensitive, and you can use any of the advanced features detailed in the advanced page of [BetterInputFiles.jl](https://www.omegalambda.au/BetterInputFiles.jl/dev/advanced/). The following is an example of creating a covariance matrix for the DES, SNLS, and SDSS instruments.
```toml
[ CovarianceMatrix ]

    [ CovarianceMatrix.Instrument.DES ]
        Astrometry = 0.001
        NonLinearity = 0.001
        PhotometricZeroPointing = 0.002
        PhotometricBias = 0.003
        Uniformity = 0.00381 # 0.0066 / √3
        Filters =               ["g",   "r",    "i",    "z"]
        ABUncertainty =         [0.002, 0.002,  0.001,  0.005]
        SuperCalCorrection =    [0.0,   0.0,    0.0,    0.0]
        FilterUncertainty =     [0.6,   0.6,    0.6,    0.6]
        CentralWavelength =     [481,   644,    781,    913]

    [ CovarianceMatrix.Instrument.SNLS ]
        Filters =               ["g",   "r",    "i",    "z"]
        ZPUncertainty =         [0.002, 0.002,  0.002,  0.002]
        SuperCalCorrection =    [0.007, -0.001, -0.006, 0.002] 
        FilterUncertainty =     [0.3,   1.0,    3.1,    0.6]
        CentralWavelength =     [475,   640,    766,    925]

    [ CovarianceMatrix.Instrument.SDSS ]
        Filters =               ["u",   "g",    "r",    "i",    "z"]
        ZPUncertainty =         [0.023, 0.002,  0.002,  0.002,  0.002]
        SuperCalCorrection =    [0.0,   -0.003, 0.004,  0.001,  -0.008]
        FilterUncertainty =     [0.7,   0.6,    0.6,    0.6,    0.6]
        CentralWavelength =     [355,   487,    617,    748,    893]
```
`FilterUncertainty` and `CentralWavelength` are assumed to be in `nm`, and `ABUncertainty`, `ZPUncertainty`, and `SuperCalCorrection` are assumed to be in `mag`. This will produce a `CovarianceMatrix` object, and save it in a `.jld2` file.

## Analysis
You can also do some basic analysis, such as plotting and drawing random offsets from the covariance matrix. The latter is particularly useful when combined with [SALTJacobian.jl](https://github.com/dessn/SALTJacobian.jl.git), allowing the creation of arbitrary, covariantly random SALT surfaces. 
```toml
[ CovarianceMatrix ]
    input = "../Outputs/DES/DES.jld2" # Import premade covariance matrix

[[ Analaysis ]]
    output = "Plots" # Will create a new directory inside global / output_path
    [ Analysis.plot ]
        # Plot options go here

[[ Analysis ]]
    output = "Draw"
    [ Analysis.draw ]
        num = 100 # Draw 100 random offsets
        template = "saltjacobian.toml" # Place the 100 random offsets into this SALTJacobian template
```
