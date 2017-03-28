"""
Computes the residual of asset market for the histogram approximation of distribution.
```
    K′ - ∫ a′() dΦ
```
The law of motion for histogram is done as in Young (2010);

### INPUT
- `cp::ConsumerProblem`
- `K` : capital guess
- `ss::StstHistogram`

### RETURNS
* residual if ss is `Void`
* update `ss` otherwise
"""
function stst_histogram_resid(cp::ConsumerProblem, K::Real, ss::Union{StstHistogram,Void} = Void())

    nEps, nSavingsPar = cp.nEps, cp.nSavingsPar

    #== Recover prices ==#
    R      = 1 + netintr(K);
    wage   = wagefunc(K);

    #=======================================================#
    ### Compute Policy function for set of prices         ###
    #=======================================================#
    @printf(" Inner loop (Policy) \n")
    f!(Θ, fvec) = eulerres2!(fvec, Θ, Θ, R, R, wage, wage, cp)
    res = nlsolve(f!, cp.Θinit, autodiff = true)
    @printf("   - Policy converged to %.2e \n", res.residual_norm)

    #== Save and reshape solution ==#
    copy!(cp.Θinit, res.zero)
    mΘ = reshape(res.zero, nSavingsPar, nEps)

    #=================================================================#
    ### Compute Stationary Distribution from decision rules         ###
    #=================================================================#

    #== Compute Transition matrix ==#
    Πaggr = forward_mat(cp, mΘ)

    #== Invariant distribution ==#
    _, eigenv = eigs(Πaggr, nev = 1,which=:LM);
    eigenv = real( squeeze(eigenv,2) )
    vHistogram = eigenv ./ sum(eigenv);

    #== CHECK all elements are positive ==#
    minimum(vHistogram) > - 1e-12 || throw(error("Negative histogram weight"))
    vHistogram[vHistogram .< 1e-12] = 0.0
    vHistogram = vHistogram ./ sum(vHistogram)

    #== CHECK consistency ==#
    err_invDist = maximum( abs( Πaggr*vHistogram - vHistogram) );
    @printf("   - Invariant distribution Converged to %.2e \n", err_invDist)
    err_invDist<1e-8 || throw(error("Invariante distribution is incorrect"))

    #== Expected Capital ==#
    Ksupply = expect_k(vHistogram, cp)
    @printf("   * K       : %6.4f\n", K)
    @printf("   * Ksupply : %6.4f\n", Ksupply)
    @printf("   * resid   : %6.4f\n", Ksupply-K)
    println("---------------------------------------")

    if isa(ss, Void)
        return K - Ksupply
    else
        #== Update ss variable ==#``
        copy!(ss.mΘ, mΘ)
        copy!(ss.vHistogram, vHistogram)
        copy!(ss.mHistogram, vHistogram)

        ss.R     = R
        ss.wage  = wage
        ss.capital = Ksupply


        xstst = [vHistogram[2:end]; Ksupply; 0.0] ;
        ystst = vec(mΘ);
        copy!(ss.xstst, xstst)     # [xHistogram, K, ϵ]
        copy!(ss.ystst, ystst)                                    # [Θ]
        copy!(ss.Zstst, [xstst; ystst; xstst; ystst; 0.0])

        return Void
    end
end

# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

"""
Computes residual of asset market using the density as an approximation
for the distribution.

### INPUT
    - `cp::ConsumerProblem`
    - `K`
    - `mMoments`
    - `init_moments`
    - `ss_density`
"""
function stst_density_resid(cp::ConsumerProblem, K::Float64, init_moments::Matrix{Float64}, ss_density::Union{Void,StstDensity} = Void())

    nSavingsPar = cp.nSavingsPar
    nEps        = cp.nEps

    #== Recover prices ==#
    R      = 1 + netintr(K);
    wage   = wagefunc(K);

    #=======================================================#
    ###  Compute Policy function for set of prices        ###
    #=======================================================#
    @printf("  Inner loop (Policy) \n")
    f!(Θ, fvec) = eulerres2!(fvec, Θ, Θ, R, R, wage, wage, cp)
    res = nlsolve(f!, cp.Θinit, autodiff = true)
    @printf("   - Policy converged to %.2e \n", res.residual_norm)
    #== Save and reshape solution ==#
    copy!(cp.Θinit, res.zero)
    mΘ = reshape(res.zero, nSavingsPar, nEps)

    #==============================================================#
    ###  Compute Stationary Distribution from decision rules     ###
    #==============================================================#
    nMoments = size(init_moments,1)
    nCoeff   = nMoments+1
    nQuad    = cp.nQuad

    asset_quad_points, quad_weights = cp.asset_quad_points, cp.quad_weights

    #== COMPUTE policy on quad nodes ==#
    mAssetPrimeonQuad = zeros(nQuad, nEps)
    for ieps=1:nEps

        S = interpolate(cp, mΘ[:,ieps])
        mAssetPrimeonQuad[:,ieps] = S[asset_quad_points]
    end

    #== Find parameters from Initial moments ==#
    mDensityCoeff  = zeros(nCoeff, nEps)
    mParamCoeff    = zeros(nMoments, nEps)

    mDensity = zeros(nQuad,nEps)

    mMoments       = copy(init_moments)
    mMomentsNext   = similar(init_moments)
    mMomentsCheck  = similar(init_moments)

    options = Optim.Options(g_tol = 1e-12, iterations = 100)

    iter = 1
    toldensity = 1e-5
    err = 1.0
    max_iterdens = 2000
    while (err>toldensity) & (iter<=max_iterdens)

        #------------------------------------------------------#
        ##  FIND coefficients of density GIVEN moments        ##
        #------------------------------------------------------#
        for ieps = 1:nEps

            #== Set up min PROBLEM ==#
            f(ρ) = density_int(ρ, mMoments[:,ieps], cp)
            g!(ρ,∇dens) = density_int_g!(ρ, ∇dens, mMoments[:,ieps], cp)
            res = optimize(f, g!, mParamCoeff[:,ieps])

            #== Solution ==#
            ρ_star        = res.minimizer;
            normalization = res.minimum;
            mParamCoeff[:,ieps] = ρ_star;

            # storage parameters
            mDensityCoeff[:,ieps] = [ 1. / normalization; ρ_star]

            # storage density values
            mDensity[:,ieps] = 1 / normalization * density(ρ_star, mMoments[:,ieps], asset_quad_points, cp)

            for iMom = 1 : nMoments
                mMomentsCheck[iMom,ieps] = density_int(ρ_star, mMoments[:,ieps], cp, iMom, 1 / normalization )
            end

        end

        #== CONSITENCY CHECK ==#
        err_moments = maximum( abs( mMoments[:] - mMomentsCheck[:] ));
        # @printf("    Iteration %3d:  Moments t CHECK: %.3e ",iter, err_moments)

        #-------------------------------------------#
        ##  NEXT period implied Moments            ##
        #-------------------------------------------#
        for ieps = 1 : nEps # WARNING NEXT period shock

            ##  NOTE:  computing next period first moments     ##
            ## m' = ∑ π( ̃ϵ | ϵ ) ∫ a'( ̃ϵ , a; z, m ) g_ϵ(a) da
            law_of_motion_aux = 0.0;

            for jeps = 1 : nEps # WARNING THIS period shock

                law_of_motion_aux = law_of_motion_aux +
                        cp.z_dist[jeps] * cp.Π[jeps, ieps] *
                         dot( quad_weights,  mAssetPrimeonQuad[:,jeps] .* mDensity[:, jeps] )
            end
            mMomentsNext[1,ieps] = law_of_motion_aux / cp.z_dist[ieps];

            # Compute higher order moments (centered)
            for iMom = 2 : nMoments

                ##  NOTE:  computing next period ith moment                    ##
                ##         Moment inside is with respect to ieps               ##
                ##
                ##     mᵢ'(ϵ) = ∑ π( ̃ϵ | ϵ ) ∫ [ a'( ̃ϵ , a; z, m ) - m¹'(ϵ) ]^i g_̃ϵ (a) da
                ##
                law_of_motion_aux = 0.0

                for jeps = 1 : nEps # WARNING THIS period shock

                    law_of_motion_aux = law_of_motion_aux +
                        cp.z_dist[jeps] * cp.Π[jeps, ieps] *
                        dot(quad_weights,
                         ( ( mAssetPrimeonQuad[:, jeps] - mMomentsNext[1,ieps] ) .^ iMom ) .* mDensity[:, jeps] );
                end

                mMomentsNext[iMom, ieps] = law_of_motion_aux / cp.z_dist[ieps];
            end

        end

        #== Distance between moments ==#
        err = maximum( abs( mMomentsNext[:] - mMoments[:] ) );
        # print iterations.
        # @printf(" distance: %.3e \n", err)

        #== UPDATE ==#
        mMoments[:]       = mMomentsNext[:]

        if (err<=toldensity) || (max_iterdens==500)
            @printf("   - Solving Invariant density:  Iter   Mom t CHECK   distances\n")
            @printf("                                 %3d     %.3e     %.3e \n", iter, err_moments, err)

        else

            mDensityCoeff[:]  = 0.0
            mParamCoeff[:]    = 0.0
            mDensity[:]       = 0.0
            mMomentsNext[:]   = 0.0
            mMomentsCheck[:]  = 0.0
        end

        iter += 1
    end
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #== Expected Capital ==#
    Ksupply = dot(cp.z_dist[:], mMoments'[:,1])
    @printf("   * K       : %6.4f\n", K)
    @printf("   * Ksupply : %6.4f\n", Ksupply)
    @printf("   * resid   : %6.4f\n", Ksupply-K)
    println("---------------------------------------")

    if isa(ss_density, Void)
        return K - Ksupply
    else
        #== Update ss_density variable ==#``
        copy!(ss_density.mΘ, mΘ)
        copy!(ss_density.mMoments, mMoments)
        copy!(ss_density.mDensityCoeff, mDensityCoeff)

        ss_density.R     = R
        ss_density.wage  = wage
        ss_density.capital = Ksupply

        return Void
    end
end

# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

"""
Function to compute expectation
"""
function expect_k{T<:Real}(vHistogram::Vector{T}, cp::ConsumerProblem)

    asset_grid_fine = cp.asset_grid_fine
    nEps = size(cp.Π,1); nHistogram = length(vHistogram)

    vHistogramK = reshape(vHistogram, div(nHistogram,nEps), nEps)
    vHistogramK = sum(vHistogramK,2)

    return dot(squeeze(vHistogramK,2), asset_grid_fine)
end

"""
Compute the Transition matrix from policies

"""
function forward_mat{T<:Real}(cp::ConsumerProblem, Θ::Array{T})

    Π, asset_grid_fine = cp.Π, cp.asset_grid_fine

    #== Sizes ==#
    nAssetsFine = cp.nAssetsFine
    nSavingsPar, nEps = cp.nSavingsPar, cp.nEps

    #== Reshape policies ==#
    mΘ = reshape(Θ, nSavingsPar, nEps)

    #== Initialize Vectors for sparse matrix ==#
    Ic  = Array(Int64,0)
    Ir  = Array(Int64,0)
    Val = Array(T,0)         # should contain Reals

    for ieps = 1:nEps

        #== Saving policy at t ==#
        Sthis::Interpolations.Extrapolation{T,1} = interpolate(cp, mΘ[:,ieps])

        #== Assets today ==#
        Athis = asset_grid_fine

        #== Policy at grid ==#
        Anext::Vector{T} = Sthis[Athis]

        ic, ir, vv = linear_trans(Athis, Anext)

        for jeps = 1:nEps
            ptrans = Π[ieps, jeps]
            ioff = (ieps-1) * nAssetsFine
            joff = (jeps-1) * nAssetsFine

            #== Append on indices ==#
            append!(Ic, ic + ioff)
            append!(Ir, ir + joff)
            append!(Val, ptrans * vv )
        end
    end
    nΠaggr = nAssetsFine * nEps

    return sparse(Ir, Ic, Val, nΠaggr, nΠaggr)
    # return Ir, Ic, Val
end

"""
Computes the ...

### OUTPUT
- `iFr`  :
- `iTo`  :
- `prob` :

### Comments
- ``
"""
function linear_trans{T<:Real}(grid::Vector{Float64}, pol::Vector{T})

    #== Initialize vectors ==#
    iFr   = collect( 1:length(grid) )
    iTo   = similar(iFr)
    pHigh = Array(T, size(iFr))                                  ## IMPORTANT : able to hold Dual

    for (i,value) in enumerate(pol)
        # println("($i,$value)")
        ipos   = max( searchsortedfirst(grid, value), 2)         ## NOTE: 1st value ≥ in grid

        iTo[i]   = ipos
        pHigh[i] = ( value - grid[ipos-1] ) / ( grid[ipos] - grid[ipos-1] )
    end

    #== Repeat indices ==#
    iFr  = [iFr; iFr ]
    iTo  = [iTo; iTo-1 ]
    prob = [pHigh; 1.0 - pHigh]

    return iFr, iTo, prob
end
