
"""
Holds all equilibrium conditions. Used to do the linearization Step

### INPUT
- `Y`       : [X_t, X_{t-1}, η_t, ϵ_{t-1}]
              Variable X at t, t-1, expectational shocks ans regular shocks
              where X = [vHistogramDev; Kaggr; dx; vSavingsPar]
- `iVarY`   : index

"""
function equil_histogram{T<:Real}(Y::Vector{T}, iVarY::Dict, ss::StstHistogram, cp::ConsumerProblem)

    #== Extract variables ==#
    X     = Y[ iVarY[:X]    ]
    Xlag  = Y[ iVarY[:Xlag] ]
    eta   = Y[ iVarY[:eta]  ]
    eps   = Y[ iVarY[:eps]  ]

    #== Create fnc to unpack ==#
    vHistogram, Θ, dz, KAggr             = unpack_histogram(X,ss)
    vHistogramlag, Θlag, dzlag, KAggrlag = unpack_histogram(Xlag,ss)

    #== Recover PRICES ==#
    R      = 1 + netintr(KAggr, dz);
    wage   = wagefunc(KAggr, dz);

    Rlag      = 1 + netintr(KAggrlag, dzlag);
    wagelag   = wagefunc(KAggrlag, dzlag);

    # ******************************************************************
    #   EQUATIONS
    # ==================================================================


    # Dynamics of distribution of wealth
    Πaggr = forward_mat(cp, Θlag)

    vHistogram_lom = Πaggr * vHistogramlag

    resD = distr2x(vHistogram_lom, ss.vHistogram) - distr2x(vHistogram, ss.vHistogram)

    # Exogenoous Shocks
    resZ = dz - ρz * dzlag - σz * eps

    # Aggregate Capital
    resK = KAggr - expect_k(vHistogram, cp);

    # Household policy rules
    resC = eulerres( Θlag, Θ, Rlag, R, wagelag, wage, cp) + eta


    resid = [resD;          # Distribution
             resZ;          # Exogenous shocks
             resK;          # Aggregate capital
             resC]          # Household policies

end

# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

"""
Unpack the vector X (variable values) into its
components (distribution of wealth, policy rule parameters, aggregate
variables)
"""
function unpack_histogram{T<:Real}(X::Vector{T}, ss::StstHistogram)

    ## HISTOGRAM ##
    vHistogram = x2distr( X[ss.iXvar[:histogram]], ss.vHistogram );

    ## SHOCK ##
    dz = X[ss.iXvar[:aggr_exo]];

    ## AGGREGATE VAR ##
    KAggr = X[ss.iXvar[:aggr_end]];

    ## SavingsPar ##
    Θ = X[ss.iXvar[:household]];

    return vHistogram, Θ, dz, KAggr
end

function x2distr{T<:Real}(vHistogramDev::Vector{T}, vHistogramSs::Vector{Float64})

    vHistogram = Array(T, length(vHistogramSs) )
    vHistogram[:]  = vHistogramSs + [-sum(vHistogramDev); vHistogramDev]
end

function x2distr!{T<:Real}(vHistogram::Vector{T}, vHistogramDev::Vector{T}, vHistogramSs::Vector{Float64})

    copy!(vHistogram,vHistogramSs + [-sum(vHistogramDev); vHistogramDev])
end

function distr2x{T<:Real}(vHistogram::Vector{T}, vHistogramSs::Vector{Float64})

    vHistogramDev = vHistogram - vHistogramSs

    return vHistogramDev[2:end]
end
