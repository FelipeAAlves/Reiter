
"""
Holds all equilibrium conditions. Used to do the linearization Step

### INPUT
- `Y`       : [X_t, X_{t-1}, η_t, ϵ_{t}]
              Variable X at t, t-1, expectational shocks ans regular shocks
              where X = [vHistogramDev; capital; dx; vSavingsPar]
- `iVarY`   : index

"""
function equil_histogram{T<:Real}(Y::Vector{T}, iZvar::Dict, ss::StstHistogram, cp::ConsumerProblem)

    #== Extract variables ==#
    x′::Vector{T}   = Y[ iZvar[:x′] ]
    y′::Vector{T}   = Y[ iZvar[:y′] ]
    x::Vector{T}    = Y[ iZvar[:x] ]
    y::Vector{T}    = Y[ iZvar[:y] ]
    ϵ_shocks::Vector{T}    = Y[ iZvar[:eps]  ]

    #== Create fnc to unpack ==#
    vHistogram′, K′, dZ′  = unpack_x(x′,ss.ixvar)
    vHistogram , K , dZ   = unpack_x(x ,ss.ixvar)
    Θ′ = y′;
    Θ  = y ;
    #== Recover PRICES ==#
    R      = 1 + netintr(K, dZ);
    wage   = wagefunc(K, dZ);

    R′      = 1 + netintr(K′ , dZ′ - σz * ϵ_shocks[1]);
    wage′   =     wagefunc(K′, dZ′ - σz * ϵ_shocks[1]);

    # ******************************************************************
    #   EQUATIONS
    # ==================================================================

    ## Equation Household policy rules  ##
    resC = eulerres( Θ, Θ′, R, R′, wage, wage′, cp)

    ## Equation Dynamics of distribution of wealth ##
    Πaggr = forward_mat(cp, Θ)

    vHistogram_lom = Πaggr * vHistogram

    resD = distr2x(vHistogram′) - distr2x(vHistogram_lom)

    ## Equation Exogenoous Shocks ##
    resZ = dZ′ - ρz * dZ - σz * ϵ_shocks[1]

    ## Equation Aggregate Capital  - WARN make it end of period capital (otherwise, problem with the state) ##
    resK = K′ - expect_k(vHistogram′, cp);



    resid = [resC;          # Household policies
             resD;          # Distribution
             resK;          # Aggregate capital
             resZ]          # Exogenous shocks

end

# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

"""
Unpack the vector x (variable values) into its
components (distribution of wealth, aggregate
variables)
"""
function unpack_x{T<:Real}(X::Vector{T}, ixvar::Dict{Symbol,UnitRange{Int64}})

    ## HISTOGRAM ##
    vHistogram = x2distr( X[ ixvar[:histogram] ]);

    ## AGGREGATE Endo VAR ##
    capital = X[ixvar[:aggr_end]][1];

    ## SHOCK ##
    dZ = X[ixvar[:aggr_exo]][1];

    return vHistogram, capital, dZ
end

function x2distr{T<:Real}(xhistogram::Vector{T})

    vHistogram = Array(T, length(xhistogram)+1 )
    copy!(vHistogram, [1.0-sum(xhistogram); xhistogram])
end

function distr2x{T<:Real}(vHistogram::Vector{T})

    return vHistogram[2:end]
end
