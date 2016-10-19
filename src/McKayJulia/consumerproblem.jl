"""
Description of the Household Problem

## Fields
- `β`
- `γ`

- `a̲`
- `euler_asset_knots`
- `asset_grid_fine`

- `z_vals`
- `Π`
"""
type MckayHousehold

    ##  Parameters  ##
    β::Float64
    γ::Float64

    ##  Asset Grid (Policies) ##
    a̲::Float64

    nSavingsPar::Int64
    euler_asset_knots::Vector{Float64}

    nLaborPar::Int64
    labor_asset_grid::Vector{Float64}

    ##  Asset Grid (Distribution) ##
    asset_grid_fine::Vector{Float64}

    ##  Quadrature  ##
    quad_weights::Vector{Float64}
    asset_quad_points::Vector{Float64}
    nQuad::Int64
    nMoments::Int64

    ##  Exogenous Shocks  ##
    z_vals::Vector{Integer}
    Π::Matrix{Float64}
    z_dist::Vector{Float64}
    nEps::Int64

    ##  Fiscal Policy  ##
    μ::Float64
    τ::Float64

    ##  EXTRA NOTE  ##
    indSavingsPar::UnitRange{Int64}
    indLaborPar::UnitRange{Int64}
    Θinit::Vector{Float64}
    KrepSS::Float64
end

"""
Constructor of MckayHousehold
"""
function MckayHousehold()

    #== Utility Parameters ==#
    β = 0.96
    γ = 1.0

    KrepSS = ( ( α * (aggEmployment ^ (1 - α)) ) / ( (1 / β) - (1 - δ) ) ) ^ (1 / (1 - α))

    ##  DISCRETIZATION of Income Process  ##
    # *************************************************************************************
    nEps = 2                                                                    # permanent productivity states
    z_vals = [0 ; 1]
    uDuration = 1                                                               # number of period unemployed before employed

    #== Transition Probabilities ==#
    Π = zeros(nEps,nEps)
    Π[1,:] = [uDuration/(1 + uDuration) 1-(uDuration / (1 + uDuration))]
    Π[2,1] = ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration)))
    Π[2,2] = 1-((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration)))

    z_dist = [1 - aggEmployment; aggEmployment];                         # distribution

    μ = .15;                                                                    # unemployment benefits
    τ = μ * (1.0 - aggEmployment) / aggEmployment;                              # flat tax
    # -------------------------------------------------------------------------------------

    # ******************************************************************
    #   IMPLEMENTATION CHOICES
    # ==================================================================

    # LIMITS on Capital capital:
    a̲ = 0.0
    assets_min = a̲;
    assets_max = 3 * KrepSS;                               # original = 3, CHANGE line 127 if increased to have 1.5 KrepSS
                                                           #               CHANGE nAssetsQuad to compute integral
    #AssetsMin = aabar;  AssetsMax = 5;

    # Number of nodes in Policies
    nSavingsPar    = 120;                                                       # policy function; default = 80-160
    indSavingsPar  = 1:nSavingsPar

    nLaborPar      = 100;
    indLaborPar    = nSavingsPar+1:nSavingsPar+nLaborPar

    nHistogram     = 75;                        # finer grid     ; default = 75-100
                                                #       bigger values perform even worst during linearization
                                                #       bc what changes is only the weigth bet adjancent points but not the POINT where I go
    # APPROXIMATION of density
    nMoments     = 5;                                                        # original: 3
    nQuad        = 12;                                                       # original: 8,
    #                                                                             #           12 works for AssetsMax = 3 * KRepSS
    #                                                                             #           16 works for AssetsMax = 4 * KRepSS
    #                                                                             #           20 works for AssetsMax = 5 * KRepSS
    # nMomentsTotal = nEps * nMoments;                                           # total moments
    # nDensityCoeffTotal = nEps * nDensityCoeff;                                 # total coeff

    # ******************************************************************
    #   FINE grid for wealth histogram:
    # ==================================================================
    nHistogramTotal = nHistogram * nEps   # makes only sense for nEps==1; otherwise_:DSF!
    #
    # Accumulate points in the beginning
    # [AssetsGridFine, logshift] = makeknotd(AssetsMin,AssetsMax,nHistogram);
    # [temp, logshift] = makeknotd(AssetsMin,AssetsMax,nHistogram);
    #
    # Equidistant GRID
    asset_grid_fine = collect( linspace(assets_min, assets_max, nHistogram) )

    # ******************************************************************
    #   Quadature weights
    # ==================================================================
    quad_points, quad_weights = gausslegendre( nQuad )
    asset_quad_points = scale_up(quad_points, assets_min+1e-1, assets_max)

    # ******************************************************************
    #   knot points for EULER Collocation
    # ==================================================================
    xmin = 0.000
    xmax = assets_max * 3/8 # beyond that point, linear extrapolation is fine!

    n1        = ceil( nSavingsPar/3 )
    x1        = linspace(xmin, 0.25, n1+1)
    x2        = logspaceshift(0.25, xmax, nSavingsPar-n1, 1.0)

    euler_grid  = [collect(x1[1:end-1]); x2]
    euler_asset_knots = euler_grid[2:end]                                             # take out the zero

    # ******************************************************************
    #   Knot Points for LABOR SUPPLY POLYNOMIAL:
    # ==================================================================
    nxmin = 0.001;
    nxmax = assets_max;
    n1 = ceil(nLaborPar/3);
    x1 = linspace(nxmin,0.2,n1+1);
    x2 = logspaceshift(0.2, nxmax, nLaborPar-n1, 1.0);

    labor_asset_grid = [collect(x1[1:end-1]); x2];
    # -------------------------------------------------------------------------------------


    # ******************************************************************
    #   Steady-State NOTErmation
    # ==================================================================

    #== Initial GUESS for the policy coefficients ==#
    Par0 = [0.0; 0.80*euler_asset_knots];                                             # 0.8 seems to work
    Θinit = repmat(Par0, length(z_vals))


    MckayHousehold(β, γ, a̲,
                    nSavingsPar, euler_asset_knots,
                    nLaborPar, labor_asset_grid,
                    asset_grid_fine,
                    quad_weights, asset_quad_points, nQuad, nMoments,
                    z_vals, Π, z_dist, nEps, μ, τ,
                    indSavingsPar, indLaborPar, Θinit, KrepSS)
end

# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

## TO DO : is this the best alternative to hold the Stst
"""
Type to hold NOTErmation on the Steady-State

## Fields
- `mΘ`
- `vHistogram` : histogram with weights

- `R`
- `wage`
"""
type StstHistogram
    mΘ::Matrix{Float64}
    vHistogram::Vector{Float64}
    mHistogram::Matrix{Float64}

    R::Float64
    wage::Float64

    Kaggr::Float64

    ##  Dictionary  ##
    iXvar::Dict{Symbol,Union{UnitRange{Int64}, Int64}}
end
"""
Constructor for type StstHistogram
"""
function StstHistogram(cp::MckayHousehold)

    #== Sizes ==#
    nSavingsPar, nEps = div(length(cp.Θinit),size(cp.Π,1)), size(cp.Π,1)
    nAssetsFine   = length(cp.asset_grid_fine)
    nHistogram    = nAssetsFine*nEps

    mΘ = zeros(nSavingsPar, nEps)
    vHistogram = zeros(nHistogram)
    mHistogram = zeros(nAssetsFine,nEps)

    ###  Create the Dicitionary  ###
    # *************************************************************************************
    iXvar = Dict{Symbol,Union{UnitRange{Int64}, Int64}}()

    nn = 0
    iXvar[:histogram] = 1:(nHistogram-1);                      nn = (nHistogram-1)
    iXvar[:aggr_exo]  = nn+1;
    iXvar[:aggr_end]  = nn+2;                                  nn = nn+2;
    iXvar[:household] = nn+1:(nn) + length(mΘ);
    # -------------------------------------------------------------------------------------

    StstHistogram(mΘ, vHistogram, mHistogram, 1.0, 0.0, 0.0, iXvar)
end

# **********************************************************************************************************************
#
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
#
# **********************************************************************************************************************

"""
Utility function
"""
function utilfn(cp::MckayHousehold, cons)
    cp.γ==1 ? log(cons) : cons^(1-cp.γ)
end

"""
Marginal utility function
"""
function dutil(cp::MckayHousehold, cons)

    if any( cons .< 0.0 )
        throw( throw(ArgumentError("Negative Consumption")) )
    end

    cons .^ (-cp.γ)

end

"""
Inv Marginal utility function
"""
function inv_dutil(cp::MckayHousehold, dutil)

    if any( dutil .< 0.0 )
        throw( throw(ArgumentError("Negative marginal utility of cons")) )
    end

    dutil .^ (-1.0/cp.γ)

end


"""
Description:
Grid of n points between xa and xb, equidistant in log(x+xshift)
"""
function logspaceshift(xa,xb,n,xshift);

  grid = exp( linspace(log(xa+xshift),log(xb+xshift),n) ) - xshift

  grid[1]   = xa
  grid[end] = xb

  return grid
end
