module Reiter

using Interpolations
using ForwardDiff
using NLsolve: nlsolve
using Optim
using FastGaussQuadrature

export
# Constructor
    ConsumerProblem, StstHistogram,

# xxxx
    netintr, wagefunc,
    eulerres!, eulerres_matlab!, stst_histogram_resid,
# xxxx
    test_jaco_interp!

# xxxx


##### includes
include("KrusselSmithJulia\\init_model_params.jl")
include("KrusselSmithJulia\\consumerproblem.jl")
include("KrusselSmithJulia\\eulerres.jl")
include("KrusselSmithJulia\\param_density.jl")
include("KrusselSmithJulia\\steadystate.jl")
include("KrusselSmithJulia\\equil_histogram.jl")

end # module
