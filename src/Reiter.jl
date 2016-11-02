module Reiter

import Base: show
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
include("KrusellSmithJulia/init_model_params.jl")
include("KrusellSmithJulia/consumerproblem.jl")
include("KrusellSmithJulia/eulerres.jl")
include("KrusellSmithJulia/param_density.jl")
include("KrusellSmithJulia/steadystate.jl")
include("KrusellSmithJulia/equil_histogram.jl")

### SDPricing
# using CompEcon
# using QuantEcon
# using NLsolve: nlsolve
#
# export
# # Constructor
#     FirmProblem, CollocStruct
# # xxxx
#
# ##### includes
# include("SDPricing\\firm_problem.jl")
# include("SDPricing\\solve_col.jl")
# # include("SDPricing/firm_problem.jl")
# # include("SDPricing/solve_col.jl")

end # module
