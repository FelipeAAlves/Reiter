
"""
Holds information of Household Problem

### Fields

- `β` : discount factor
- `σ` : risk aversion

- `a̲`           : borrowing constraint
- `nSavingsPar` : policy discretization
- `asset_knots` : policy discretization

- `nAssetsFine`     : histogram points
- `asset_grid_fine` : histogram points

- `nQuad`   : quadrature info (used in density only)
- `nMoments`: quadrature info (used in density only)

- `z_vals` : exogenous shocks
- `Π`      : transition matrix
"""
type ConsumerProblem

    ##  Parameters  ##
    β::Float64
    σ::Float64

    ##  Asset Grid (Policies) ##
    nSavingsPar::Int64
    a̲::Float64
    asset_knots::Vector{Float64}

    ##  Asset Grid (Distribution) ##
    nAssetsFine::Int64
    asset_grid_fine::Vector{Float64}

    ##  Quadrature  ##
    nQuad       ::Int64
    nMoments    ::Int64
    quad_weights::Vector{Float64}
    asset_quad_points::Vector{Float64}

    ##  Exogenous Shocks  ##
    nEps::Int64
    z_vals::Vector{Integer}
    Π::Matrix{Float64}
    z_dist::Vector{Float64}

    ##  Fiscal Policy  ##
    μ::Float64
    τ::Float64

    ##  Extra NOTE  ##
    Θinit::Vector{Float64}
    KrepSS::Float64
end

"""
Constructor of ConsumerProblem
"""
function ConsumerProblem()

    #== Utility Parameters ==#
    β = 0.96
    σ = 1.0

    KrepSS = ( ( α * (aggEmployment ^ (1 - α)) ) / ( (1 / β) - (1 - δ) ) ) ^ (1 / (1 - α))

    # ******************************************************************
    #   DISCRETIZATION of Income Process
    # ==================================================================
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

    # number of grid points for each shock:
    nSavingsPar    = 120;                       # policy function; default = 80-160
    nAssetsFine    = 200;                       # finer grid     ; default = 75-100
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
    nHistogram = nAssetsFine * nEps   # makes only sense for nEps==1; otherwise_:DSF!
    #
    # Accumulate points in the beginning
    # [AssetsGridFine, logshift] = makeknotd(AssetsMin,AssetsMax,nAssetsFine);
    # [temp, logshift] = makeknotd(AssetsMin,AssetsMax,nAssetsFine);
    #
    # Equidistant GRID
    asset_grid_fine = collect( linspace(assets_min, assets_max, nAssetsFine) )

    # ******************************************************************
    #   Quadature weights
    # ==================================================================
    quad_points, quad_weights = gausslegendre( nQuad )
    asset_quad_points = scale_up(quad_points, assets_min+1e-3, assets_max)
    quad_weights = 0.5 * (assets_max-(assets_min+1e-3)) * quad_weights

    # ******************************************************************
    #   knot points for EULER Collocation
    # ==================================================================
    xmin = 0.000
    xmax = assets_max * 1/2     # default 3/8 beyond that point, linear extrapolation is fine!

    n1        = ceil( nSavingsPar/3 )
    x1        = linspace(xmin, 0.25, n1+1)
    x2        = logspaceshift(0.25, xmax, nSavingsPar-n1, 1.0)

    euler_grid  = [collect(x1[1:end-1]); x2]
    asset_knots = euler_grid[2:end]                                             # take out the zero
    # -------------------------------------------------------------------------------------

    # ******************************************************************
    #   Steady-State info
    # ==================================================================

    #== Initial GUESS for the policy coefficients ==#
    Par0 = [0.0; 0.80*asset_knots];                                             # 0.8 seems to work
    Θinit = repmat(Par0, length(z_vals))


    #== Construct ==#
    ConsumerProblem(β, σ,
                    nSavingsPar, a̲, asset_knots,
                    nAssetsFine, asset_grid_fine,
                    nQuad, nMoments, quad_weights, asset_quad_points,
                    nEps, z_vals, Π, z_dist, μ, τ,
                    Θinit, KrepSS)
end

# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

## TO DO : is this the best alternative to hold the Stst?
"""
Type to hold information regarding Steady-State for the Histogram case

## Fields
- `mΘ`
- `vHistogram` : histogram with weights

- `R`
- `wage`
"""
type StstHistogram
    ##  Policies and Distribution ##
    mΘ::Matrix{Float64}
    vHistogram::Vector{Float64}
    mHistogram::Matrix{Float64}

    ##  Aggregates  ##
    R::Float64
    wage::Float64
    capital::Float64

    ##  Dictionary  ##
    xstst::Vector{Float64}
    ystst::Vector{Float64}
    ixvar::Dict{Symbol,UnitRange{Int64}}
    iyvar::Dict{Symbol,UnitRange{Int64}}

    Zstst::Vector{Float64}
    iZvar::Dict{Symbol,UnitRange{Int64}}
end
# Observation: Let it be type to facilitate UPDATE in steady-state fnc

"""
  Constructor for type StstHistogram
"""
function StstHistogram(cp::ConsumerProblem)

    #== Sizes ==#
    nSavingsPar, nEps = cp.nSavingsPar, cp.nEps
    nAssetsFine   = cp.nAssetsFine
    nHistogram    = nAssetsFine * nEps

    mΘ = zeros(nSavingsPar, nEps)
    vHistogram = zeros(nHistogram)
    mHistogram = zeros(nAssetsFine,nEps)

    ###  Create the Dicitionary IMPORTANT ###
    # *************************************************************************************
    ixvar = Dict{Symbol,UnitRange{Int64}}()
    iyvar = Dict{Symbol,UnitRange{Int64}}()

    ##  CASE:  State variables  ## # Φ, K, z
    nx = 0
    ixvar[:histogram] = 1:(nHistogram-1);                      nx = (nHistogram-1)
    ixvar[:aggr_end]  = nx+1:nx+1;
    ixvar[:aggr_exo]  = nx+2:nx+2;
    nx += 2
    ##  CASE:  controls  ##
    iyvar[:household] = 1:length(mΘ);                         ny = length(mΘ)
    # -------------------------------------------------------------------------------------

    #=======================================================#
    ###          CONSTRUCT iZvar for use in the          ###
    #=======================================================#
    iZvar = Dict{Symbol, UnitRange{Int64}}()

    nExog = 1

    #== iZvar ==#
    iZvar[:x′]    = 1:nx
    iZvar[:y′]    = nx+1:ny+nx

    iZvar[:y] = (nx+ny) + (iZvar[:y′])
    iZvar[:x] = (nx+ny) + (iZvar[:x′])
    nn = 2*(nx+ny)
    iZvar[:eps] = nn+1:nn+nExog

    StstHistogram(mΘ, vHistogram, mHistogram, 1.0, 0.0, 0.0,
            Array(Float64,nx), Array(Float64,ny), ixvar, iyvar,
            Array(Float64,nn+nExog),iZvar)
end

"""
Type to hold information regarding Steady-State for the Histogram case

## Fields
- `mΘ`
- `vHistogram` : histogram with weights

- `R`
- `wage`
"""
type StstDensity
    ##  Policies and Distribution ##
    mΘ::Matrix{Float64}
    mMoments::Matrix{Float64}
    mDensityCoeff::Matrix{Float64}

    ##  Aggregates  ##
    R::Float64
    wage::Float64
    capital::Float64

    ##  Dictionary  ##
    iXvar::Dict{Symbol,UnitRange{Int64}}
end
"""
Constructor for type StstDensity
"""
function StstDensity(cp::ConsumerProblem)

    #== Sizes ==#
    nSavingsPar, nEps = cp.nSavingsPar, cp.nEps
    nMoments = cp.nMoments

    mΘ = zeros(nSavingsPar, nEps)
    mMoments      = zeros(nMoments, nEps)
    mDensityCoeff = zeros(nMoments+1, nEps)

    ###  Create the Dicitionary IMPORTANT ###
    # *************************************************************************************
    iXvar = Dict{Symbol,Union{UnitRange{Int64}, Int64}}()

    nn = 0
    iXvar[:density] = 1:(2*nMoments+1);                      nn = 2*nMoments+1
    iXvar[:aggr_end]  = nn+1:nn+1;
    iXvar[:aggr_exo]  = nn+2:nn+2;
    nn = nn+2;
    iXvar[:household] = nn+1:nn+length(mΘ);
    # -------------------------------------------------------------------------------------

    StstDensity(mΘ, mMoments, mDensityCoeff, 1.0, 0.0, 0.0, iXvar)
end

# **********************************************************************************************************************
#
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
#
# **********************************************************************************************************************

function scale_up(quad_points::Vector{Float64}, xmin::Float64, xmax::Float64)

    return .5 * (quad_points + 1.0) * (xmax - xmin)  + xmin
end

"""
Interest Rate function
"""
function netintr(K::Real, z::Real = 0.0)

    # adjust for differences when L ≠ 1
    KLratio = K / aggEmployment

    r = A * exp(z) * α * (KLratio) .^ (α-1) - δ

end
"""
Wage function
"""
function wagefunc(K::Real, z::Real = 0.0)

    # adjust for differences when L ≠ 1
    KLratio = K / aggEmployment;

    wage = A * exp(z) * ( 1 - α ) * (KLratio) .^ α;

end

"""
Utility function
"""
function utilfn(cp::ConsumerProblem, cons)
    cp.σ==1 ? log(cons) : cons^(1-cp.σ)
end

"""
Marginal utility function
"""
function dutil(cp::ConsumerProblem, cons)

    if any( cons .< 0.0 )
        throw( throw(ArgumentError("Negative Consumption")) )
    end

    cons .^ (-cp.σ)

end

"""
Inv Marginal utility function
"""
function inv_dutil(cp::ConsumerProblem, dutil)

    if any( dutil .< 0.0 )
        throw( throw(ArgumentError("Negative marginal utility of cons")) )
    end

    dutil .^ (-1.0/cp.σ)

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
