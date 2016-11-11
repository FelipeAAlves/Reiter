module Reiter

import Base: show
using Interpolations
using ForwardDiff
using NLsolve: nlsolve
using Optim
using FastGaussQuadrature

export
## Constructor
    ConsumerProblem, StstHistogram,
## Euler Residual
    netintr, wagefunc,
    eulerres!, eulerres_matlab!,
##
    stst_histogram_resid, stst_density_resid,
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

# TODO
#  - add an function to do expectations of Value func
#  - know now how to handle the foc - CHECK policyPerturb1dd! from Reiter
#  - introduce the Type of Variables that holds stst value and dual number structure
#  - 

# using BasisMatrices
# using QuantEcon
# using NLsolve: nlsolve
#
# export
# # Constructor
#     FirmProblem, CollocStruct
# # xxxx
#
# ##### includes
# # include("SDPricing\\firm_problem.jl")
# # include("SDPricing\\solve_col.jl")
# include("SDPricing/firm_problem.jl")
# include("SDPricing/solve_col.jl")

end # module
