module Reiter

import Base: show

# **************************************************************************************
#   Krusell-Smith ECONOMY
#
# ======================================================================================
# using Interpolations
# using ForwardDiff
# using NLsolve: nlsolve
# using Optim
# using FastGaussQuadrature

# export
# ## Constructor
#     ConsumerProblem, StstHistogram,
# ## Euler Residual
#     netintr, wagefunc,
#     eulerres!, eulerres_matlab!,
# ##
#     stst_histogram_resid, stst_density_resid,
# # xxxx
#     test_jaco_interp!
# # xxxx
#
#
##### includes
# include("tools/broydn.jl")
# include("tools/klein.jl")
# include("tools/gensys.jl")
# include("KrusellSmithJulia/init_model_params.jl")
# include("KrusellSmithJulia/consumerproblem.jl")
# include("KrusellSmithJulia/eulerres.jl")
# include("KrusellSmithJulia/param_density.jl")
# include("KrusellSmithJulia/steadystate.jl")
# include("KrusellSmithJulia/equil_histogram.jl")

# **************************************************************************************
#   SD Pricing
#
# ======================================================================================
using BasisMatrices
import BasisMatrices: nodes
using ForwardDiff
using QuantEcon
using NLsolve: nlsolve
using Roots: fzero

export
# Constructor
    StrucParameters, SetParameters

#### includes
include("tools/broydn.jl")
include("tools/klein.jl")
include("tools/gensys.jl")
include("SDPricing/pars.jl")
include("SDPricing/firm_problem.jl")
include("SDPricing/solve_col2.jl")
include("SDPricing/steadystate.jl")
include("SDPricing/equil_histogram.jl")

const __pars = SetParameters()

end # module
