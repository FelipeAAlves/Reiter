module Reiter

import Base: show
# using Interpolations
# using ForwardDiff
# using NLsolve: nlsolve
# using Optim
# using FastGaussQuadrature
#
# export
# # Constructor
#     ConsumerProblem, StstHistogram,
#
# # xxxx
#     netintr, wagefunc,
#     eulerres!, eulerres_matlab!, stst_histogram_resid,
# # xxxx
#     test_jaco_interp!
#
# # xxxx
#
#
# ##### includes
# include("KrusselSmithJulia\\init_model_params.jl")
# include("KrusselSmithJulia\\consumerproblem.jl")
# include("KrusselSmithJulia\\eulerres.jl")
# include("KrusselSmithJulia\\param_density.jl")
# include("KrusselSmithJulia\\steadystate.jl")
# include("KrusselSmithJulia\\equil_histogram.jl")

using CompEcon
using QuantEcon
using NLsolve: nlsolve

export
# Constructor
    FirmProblem, CollocStruct
# xxxx

##### includes
# include("SDPricing\\firm_problem.jl")
# include("SDPricing\\solve_col.jl")
include("SDPricing/firm_problem.jl")
include("SDPricing/solve_col.jl")

end # module
