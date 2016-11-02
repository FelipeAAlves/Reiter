"""
Computes the residual of asset market using histogram approximation of distribution.
The law of motion for histogram is done as in Young (2010);

### INPUT

"""
function stst_histogram_resid(fp:FirmProblem, w::Real, ss_histogram::Union{StstHistogram,Void} = Void())

    # nEps, nSavingsPar = cp.nEps, cp.nSavingsPar

    #=======================================================#
    ### Compute Policy function for set of prices         ###
    #=======================================================#
    @printf("   - Inner loop (Policy) \n")

    #== Save and reshape solution ==#
    Reiter.bellman_rhs!(Vnew, ξstar, coeff, w,  fp)
    
    #=================================================================#
    ### Compute Stationary Distribution from decision rules         ###
    #=================================================================#

    #== Compute Transition matrix ==#
    Πaggr = forward_mat(fp, pstar, ξstar)

    #== Invariant distribution ==#
    @time _, eigenv = eigs(Πaggr, nev = 1,which=:LM);
    eigenv = real( squeeze(eigenv,2) )
    vHistogram = eigenv ./ sum(eigenv);

    # #== CHECK all elements are positive ==#
    # minimum(vHistogram) > - 1e-12 || throw(error("Negative histogram weight"))
    # vHistogram[vHistogram .< 1e-12] = 0.0
    # vHistogram = vHistogram ./ sum(vHistogram)
    #
    # #== CHECK consistency ==#
    # err_invDist = maximum( abs( Πaggr*vHistogram - vHistogram) );
    # @printf("   - Invariant distribution Converged to %.2e \n", err_invDist)
    # err_invDist<1e-8 || throw(error("Invariante distribution is incorrect"))
    #
    # #== Expected Capital ==#
    # Ksupply = expect_k(vHistogram, cp)
    # @printf("   K       : %6.4f\n", K)
    # @printf("   Ksupply : %6.4f\n", Ksupply)
    # @printf("   resid   : %6.4f\n", Ksupply-K)
    # println("---------------------------------------")
    #
    # if isa(ss_histogram, Void)
    #     return K - Ksupply
    # else
    #     #== Update ss_histogram variable ==#``
    #     copy!(ss_histogram.mΘ, mΘ)
    #     copy!(ss_histogram.vHistogram, vHistogram)
    #     copy!(ss_histogram.mHistogram, vHistogram)
    #
    #     ss_histogram.R     = R
    #     ss_histogram.wage  = wage
    #     ss_histogram.Kaggr = Ksupply
    #
    #     return Void
    # end
end

"""
Computer the transition matrix

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
