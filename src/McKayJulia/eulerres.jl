
"""
Computes the euler and labor supply residuals at collocation nodes

### INPUTS
- `Θ`
- `Θ′`

- `R`
- `R′`

- `wage`
- `wage′`

- `π′`

- `hh::MckayHoushold`
"""
function eulerres!(out::Vector, Θ::Vector , Θ′::Vector, R::Real, R′::Real, wage::Real, wage′::Real, π′::Real, hh::MckayHoushold)

    β, Π = hh.β, hh.Π

    nEps    = hh.nEps
    nPar    = hh.nSavingsPar + hh.nLaborPar
    indS, indN    = hh.indSavingsPar, hh.indLaborPar

    #== Reshape policy coefficients ==#
    mΘ  = reshape(Θ , nPar, nEps)
    mΘ′ = reshape(Θ′, nPar, nEps)

    #== Savings and labor NEXT Period ==#
    Snext = Interpolations.AbstractInterpolation[savspline(hh, mΘ′[indS,ieps]) for ieps=1:nEps]
    Nnext = Interpolations.AbstractInterpolation[nspline(hh, mΘ′[indN,ieps]) for ieps=1:nEps]

    #== fncs to evaluate the savings & labor NEXT period ==#
    snext!(a,vals) = map!(i -> Snext[i][a], vals, 1:nEps)
    nnext!(a,vals) = map!(i -> Nnext[i][a], vals, 1:nEps)

    #== Initialize ==#
    snext = Array(Real, nEps)                                                   ##  IMPORTANT:  container must accept the Dual  ##
    nvals = Array(Real, nEps)                                                   ##  IMPORTANT:  container must accept the Dual  ##

    # loop over THIS idio shock
    for ieps = 1:nEps

        #===================================#
        ##  [1] Euler Residuasl            ##
        #===================================#

        #== Savings and Labor interpolants ==#
        Sthis, Nthis = interpolate(hh, mΘ[:,ieps])

        #== Nodes and policies decisions ==#
        athis = Sthis.itp.knots[1]; pop!(athis)                                 # pop last term (there for extrapolation ONLY)
        sthis = Sthis.itp.coefs   ; pop!(sthis)
        nthis = Nthis[athis]

        #== Recover policies ==#
        cthis = get_consumption(hh, athis, sthis, nthis, R, wage, ieps)

        #== Next period assets ==#
        anext = sthis/π′

        # loop over NEXT period assets
        for (ia, a) in enumerate(anext)

            #== Get policies NEXT period ==#
            snext!(a,snext)
            nnext!(a,nnext)

            #== Consumption NEXT period ==#
            cnext = get_consumption(hh, a, snext, nnext, R′, wage′, 1:nEps)

            if any(cnext .< 0.0) # signal to solver
                out[ia + (ieps-1)*nPar:end] = 1e2
                return Void
            else

                #== Compute expectation ==#
                MUexp = dot( dutil(hh, cnext), vec(Π[ieps,:]) )

                #== Residual ==#
                out[ia + (ieps-1)*nPar] = 1 - inv_dutil(hh, β * R′ / π′ * MUexp) / cthis[ia]
            end
        end
        # .....................................................................................

        #=====================================#
        ##  [2] Labor Residuals              ##
        #=====================================#
        athis = Nthis.itp.knots[1]; pop!(athis)
        nthis = Nthis.itp.coefs   ; pop!(nthis)
        sthis = Sthis[athis]

        #== Recover policies ==#
        cthis = get_consumption(hh, athis, sthis, nthis, R, wage, ieps)

        indUnemp = ieps .== 1 # BitArray
        labor_wage = wage * ( (1 - hh.τ) * (1 - indUnemp) )

        ind_labor_res = (indN) + (ieps-1)*nPar
        out[ind_labor_res] = ( ψ * nthis .^ ψ   ) ./ ( labor_wage * dutil(hh, cthis) ) - 1
    end

    return Void
end

function eulerres2!{T,J<:Real}(out::Vector{T}, Θ::Vector{T} , Θ′::Vector{T}, R::J, R′::J, wage::J, wage′::J, π′::J, hh::MckayHoushold)

    β, Π = hh.β, hh.Π

    nEps          = hh.nEps
    nPar          = hh.nSavingsPar + hh.nLaborPar
    indS, indN    = hh.indSavingsPar, hh.indLaborPar

    #== Reshape policy coefficients ==#
    mΘ  = reshape(Θ , nPar, nEps)
    mΘ′ = reshape(Θ′, nPar, nEps)

    #== Initialize ==#
    cthis = Array(T, nSavingsPar)                                                   ##  IMPORTANT:  container must accept the Dual  ##
    anext = Array(T, nSavingsPar)
    snext = Array(T, nSavingsPar)

    nthis = Array(T, nSavingsPar)
    nnext = Array(T, nSavingsPar)

    # loop over THIS idio shock
    for ieps = 1:nEps

        #===================================#
        ##  [1] Euler Residuasl            ##
        #===================================#

        #== Savings and Labor interpolants ==#
        Sthis, Nthis = interpolate(hh, mΘ[:,ieps])

        #== Nodes and policies decisions ==#
        athis = Sthis.itp.knots[1]::Vector{T}; pop!(athis)                                 # pop last term (there for extrapolation ONLY)
        sthis = Sthis.itp.coefs   ::Vector{T}; pop!(sthis)

        nthis = Nthis[athis]::Vector{T}

        #== Recover policies ==#
        cthis = get_consumption(hh, athis, sthis, nthis, R, wage, ieps)

        #== Next period assets ==#
        anext = sthis/π′

        # loop over NEXT period assets
        for (ia, a) in enumerate(anext)

            #== Get policies NEXT period ==#
            snext!(a,snext)
            nnext!(a,nnext)

            #== Consumption NEXT period ==#
            cnext = get_consumption(hh, a, snext, nnext, R′, wage′, 1:nEps)

            if any(cnext .< 0.0) # signal to solver
                out[ia + (ieps-1)*nPar:end] = 1e2
                return Void
            else

                #== Compute expectation ==#
                MUexp = dot( dutil(hh, cnext), vec(Π[ieps,:]) )

                #== Residual ==#
                out[ia + (ieps-1)*nPar] = 1 - inv_dutil(hh, β * R′ / π′ * MUexp) / cthis[ia]
            end
        end
        # .....................................................................................

        #=====================================#
        ##  [2] Labor Residuals              ##
        #=====================================#
        athis = Nthis.itp.knots[1]; pop!(athis)
        nthis = Nthis.itp.coefs   ; pop!(nthis)
        sthis = Sthis[athis]

        #== Recover policies ==#
        cthis = get_consumption(hh, athis, sthis, nthis, R, wage, ieps)

        indUnemp = ieps .== 1 # BitArray
        labor_wage = wage * ( (1 - hh.τ) * (1 - indUnemp) )

        ind_labor_res = (indN) + (ieps-1)*nPar
        out[ind_labor_res] = ( ψ * nthis .^ ψ   ) ./ ( labor_wage * dutil(hh, cthis) ) - 1
    end

    return Void
end


function eulerres(Θ::Vector , Θ′::Vector, R::Real, R′::Real, wage::Real, wage′::Real, hh::MckayHoushold)

    out = Array(Real, size(Θ))  ##  IMPORTANT:  container must accept the Dual  ##

    #== Compute residual ==#
    eulerres!(out, Θ , Θ′, R, R′, wage, wage′, hh)

    return out
end

# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

"""
Given a vector construct
an interpolation object that does linear interpolation in the asset dimension
"""
function Interpolations.interpolate{T<:Real}(hh::MckayHoushold, Θ::Vector{T})

    indS, indN = hh.indSavingsPar, hh.indLaborPar

    Θs = Θ[indS]; Θn = Θ[indN]

    ###                         SAVINGS Interpolant                          ###
    # **************************************************************************
    assets  = [ Θs[1]; Θs[1] + hh.asset_knots; 500.0]
    savings = [ 0.0; Θs[2:end]; 0.0 ]

    slope = ( savings[end-1] - savings[end-2] )/( assets[end-1]-assets[end-2] )
    savings[end] = savings[end-1] + ( assets[end]-assets[end-1] ) * slope

    #== Create Interpolant ==#
    S = interpolate( (assets,), savings, Gridded(Linear()) )

    #== extrapolate values outside grid ==#
    S = extrapolate(S, Flat() )
    # --------------------------------------------------------------------------

    ###                        LABOR Interpolant                             ###
    # **************************************************************************
    assets = [hh.labor_asset_grid; 500.0]
    labor  = Θn

    #== Create Interpolant ==#
    N = interpolate( (assets,), labor, Gridded(Linear()) )
    # --------------------------------------------------------------------------

    return S, N
end

function savspline{T<:Real}(hh::MckayHoushold, Θs::Vector{T})

    assets  = [ Θs[1]; Θs[1] + hh.asset_knots; 500.0]
    savings = [ 0.0; Θs[2:end]; 0.0 ]

    slope = ( savings[end-1] - savings[end-2] )/( assets[end-1]-assets[end-2] )
    savings[end] = savings[end-1] + ( assets[end]-assets[end-1] ) * slope

    #== Create Interpolant ==#
    S = interpolate( (assets,), savings, Gridded(Linear()) )

    #== extrapolate values outside grid ==#
    S = extrapolate(S, Flat() )
end

function nspline{T<:Real}(hh::MckayHoushold, Θn::Vector{T})

    assets = [hh.labor_asset_grid; 500.0]
    labor  = Θn

    #== Create Interpolant ==#
    N = interpolate( (assets,), labor, Gridded(Linear()) )
end

"""
Recover consumption
"""
function get_consumption!{T<:Real}(cons::Vector{T}, asset, savings, n, R, wage, incomeIndex::Union{AbstractVector{Int64},Integer}, hh::MckayHousehold)

    unempInd = (incomeIndex .== 1) # BitArray
    effec_wage  = wage * ( cp.μ * unempInd  + (1.0 - cp.τ) * (1 - unempInd) )

    copy!(cons, R * asset + effec_wage * n - savings)
end

function test_jaco_interp!(out::Vector, Θ::Vector, R::Real, wage::Real, hh::MckayHoushold)

    β, Π = hh.β, hh.Π

    nEps    = size(Π, 1)
    nAssets = length(hh.asset_knots) + 1

    #== Reshape policy coefficients ==#
    mΘ  = reshape(Θ, nAssets, nEps)

    # loop over idio shock
    for ieps = 1:nEps
        Sthis = interpolate(hh, mΘ[:,ieps])

        athis = Sthis.itp.knots[1]; pop!(athis)                                 # pop last term (there for extrapolation ONLY)
        sthis = Sthis.itp.coefs   ; pop!(sthis)

        out[(ieps-1)*nAssets+1:ieps*nAssets] = get_consumption(hh, athis, sthis, R, wage, ieps)
    end

end
