
"""
Holds allHolds all equilibrium conditions. Used to do the LINEARIZATION Step

### INPUT
    - `Y`  : [x', y', x,y, ϵ']
              where x = [vHistogramDev; capital; dx]
                    y = [vSavingsPar]
    - `ss` : holds some stst information
"""
function equil_histogram{T<:Real}(Y::Vector{T}, ss::StstHistogram, cp::ConsumerProblem)

    #== Extract variables ==#
    x_p1::Vector{T}   = Y[ ss.iZvar[:x′] ]
    x::Vector{T}      = Y[ ss.iZvar[:x] ]
    y_p1::Vector{T}   = Y[ ss.iZvar[:y′] ]
    y::Vector{T}    = Y[ ss.iZvar[:y] ]
    ϵ_shocks::Vector{T}    = Y[ ss.iZvar[:eps]  ]

    #== Create fnc to unpack ==#
    vHistogram_p1, K_p1, dZ_p1  = unpack_x(x_p1,ss.ixvar)
    vHistogram , K , dZ   = unpack_x(x ,ss.ixvar)
    Θ_p1 = y_p1;
    Θ  = y ;
    #== Recover PRICES ==#
    R      = 1 + netintr(K, dZ);
    wage   = wagefunc(K, dZ);

    R_p1      = 1 + netintr(K_p1 , dZ_p1 - σz * ϵ_shocks[1]);
    wage_p1   =     wagefunc(K_p1, dZ_p1 - σz * ϵ_shocks[1]);

    # ******************************************************************
    #   EQUATIONS
    # ==================================================================

    ## Equation Household policy rules  ##
    resC = eulerres( Θ, Θ_p1, R, R_p1, wage, wage_p1, cp)

    ## Equation Dynamics of distribution of wealth ##
    Πaggr = forward_mat(cp, Θ)

    vHistogram_lom = Πaggr * vHistogram

    resD = distr2x(vHistogram_p1) - distr2x(vHistogram_lom)

    ## Equation Exogenoous Shocks ##
    resZ = dZ_p1 - ρz * dZ - σz * ϵ_shocks[1]

    ## Equation Aggregate Capital  ##
    resK = K_p1 - expect_k(vHistogram_p1, cp);
    ###  WARN make it end of period capital (otherwise, problem with the state) ###


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

function x2distr!{T<:Real}(vHistogram, xhistogram::Vector{T})

    copy!(vHistogram, [1.0-sum(xhistogram); xhistogram])
end

function distr2x{T<:Real}(vHistogram::Vector{T})

    return vHistogram[2:end]
end
