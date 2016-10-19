
function eulerres{T,J <:Real}(Θ::Vector{T} , Θ′::Vector{T}, R::J, R′::J, wage::J, wage′::J, cp::ConsumerProblem)

    # out = init_matrix(Θ, Θ, Θ′, R, R′, wage, wage′)           ##  IMPORTANT:  container must accept the Dual  ##

    out = Array(T, length(Θ))
    #== Compute residual ==#
    eulerres2!(out, Θ , Θ′, R, R′, wage, wage′, cp)

    return out
end

"""
    Think if I should use this
"""
function init_matrix(x...)

    isdual = false
    for arg in x

        if isa(arg, ForwardDiff.Dual)
            isdual = true
            nder   = length(arg.partial)
        end
    end

    if isdual
        z = Array(ForwardDiff.Dual{nder,Float64}, size(x[1]) )
    else
        z = copy( x[1] )
    end
end

"""
Computes the euler residual at collocation nodes

### INPUTS
- `Θ`
- `Θ′`
- `R`
- `R′`
- `wage`
- `wage′`
- `cp::ConsumerProblem`
"""
function eulerres1!{T,J<:Real}(out::Vector{T}, Θ::Vector{T}, Θ′::Vector{T}, R::J, R′::J, wage::J, wage′::J, cp::ConsumerProblem)

    β, Π = cp.β, cp.Π

    nEps        = cp.nEps
    nSavingsPar = cp.nSavingsPar

    #== Reshape policy coefficients ==#
    mΘ  = reshape(Θ , nSavingsPar, nEps)
    mΘ′ = reshape(Θ′, nSavingsPar, nEps)

    #== Savings Next Period ==#
    Snext = Interpolations.Extrapolation{T,1}[interpolate(cp, mΘ′[:,ieps]) for ieps=1:nEps]
    snext!(a,svals) = map!(i -> Snext[i][a], svals, 1:nEps)

    #== Initialize ==#
    C     = Array(T, nSavingsPar)                                               ##  IMPORTANT:  container must accept the Dual  ##
    cnext = Array(T, nEps)
    svals = Array(T, nEps)

    # loop over idio shock
    for ieps = 1:nEps
        Sthis::Interpolations.Extrapolation{T,1} = interpolate(cp, mΘ[:,ieps])

        athis = Sthis.itp.knots[1]::Vector{T}  ; pop!(athis)                                 # pop last term (there for extrapolation ONLY)
        sthis = Sthis.itp.coefs   ::Vector{T}  ; pop!(sthis)

        get_consumption!(C, athis, sthis, R, wage, ieps, cp)

        # loop over next period assets
        anext = sthis
        for (ia, a) in enumerate(anext)

            #== Get savings next period ==#
            snext!(a,svals)
            get_consumption!(cnext, a, svals, R′, wage′, 1:nEps,cp)

            if any(cnext .< 0.0) # signal to solver
                out[ia + (ieps-1)*nSavingsPar:end] = 1e2
                return Void
            else
                #== Compute expectation ==#
                MUexp = dot( dutil(cp, cnext), vec(Π[ieps,:]) )
                cnext[:] = zero(T)

                #== Fill up Euler residual ==#
                out[ia + (ieps-1)*nSavingsPar] = 1.0 - inv_dutil(cp, β * R′ * MUexp) / C[ia]
            end
        end

        C[:] = zero(T)
    end

    return Void
end


"""
Computes the euler residual at collocation nodes (better performance)

### INPUTS
- `Θ`
- `Θ′`
- `R`
- `R′`
- `wage`
- `wage′`
- `cp::ConsumerProblem`
"""
function eulerres2!{T,J<:Real}(out::Vector{T}, Θ::Vector{T}, Θ′::Vector{T}, R::J, R′::J, wage::J, wage′::J, cp::ConsumerProblem)

    β, Π = cp.β, cp.Π

    nEps        = cp.nEps
    nSavingsPar = cp.nSavingsPar

    #== Reshape policy coefficients ==#
    mΘ  = reshape(Θ , nSavingsPar, nEps)
    mΘ′ = reshape(Θ′, nSavingsPar, nEps)

    #== Pre-allocate ==#
    cthis = Array(T, nSavingsPar)                                               ##  IMPORTANT:  container must accept the Dual  ##
    cnext = Array(T, nSavingsPar)
    anext = Array(T, nSavingsPar)
    MUexp = Array(T, nSavingsPar)
    # loop over idio shock
    for ieps = 1:nEps
        Sthis::Interpolations.Extrapolation{T,1} = interpolate(cp, mΘ[:,ieps])

        athis = Sthis.itp.knots[1]::Vector{T} ; pop!(athis)                                 # pop last term (there for extrapolation ONLY)
        sthis = Sthis.itp.coefs   ::Vector{T} ; pop!(sthis)

        get_consumption!(cthis, athis, sthis, R, wage, ieps, cp)
        copy!(anext, sthis)

        # loop over next period assets
        MUexp[:] = zero(T)
        for jeps = 1:nEps

            #== Saving policy t+1 ==#
            Snext::Interpolations.Extrapolation{T,1} = interpolate(cp, mΘ′[:,jeps]);

            #== Saving at values ==#
            snext = Snext[anext]::Vector{T}

            #== Consumption next period ==#
            get_consumption!(cnext, anext, snext, R′, wage′, jeps, cp)

            #== Marginal utility next period ==#
            MUnext = dutil(cp, cnext);
            MUexp = MUexp + Π[ieps,jeps] * MUnext;

            #== zero consumption ==#
            cnext[:] = zero(T)
        end

        #== Residual ==#
        out[ (ieps-1)*nSavingsPar+1:ieps*nSavingsPar ] = 1.0 - inv_dutil(cp, β * R′ * MUexp) ./ cthis

        #== zero consumption ==#
        cthis[:] = zero(T)
    end

    return Void
end

function eulerres_test{T,J<:Real}(Θ::Vector{T}, Θ′::Vector{T}, R::J, R′::J, wage::J, wage′::J, cp::ConsumerProblem)

    β, Π = cp.β, cp.Π

    nEps        = cp.nEps
    nSavingsPar = cp.nSavingsPar

    #== Reshape policy coefficients ==#
    mΘ  = reshape(Θ , nSavingsPar, nEps)
    mΘ′ = reshape(Θ′, nSavingsPar, nEps)

    #== Pre-allocate ==#
    cthis = Array(T, nSavingsPar)                                               ##  IMPORTANT:  container must accept the Dual  ##
    cnext = Array(T, nSavingsPar)
    anext = Array(T, nSavingsPar)
    MUexp = Array(T, nSavingsPar)
    # loop over idio shock
    # for ieps = 1:nEps
    ieps = 1
        Sthis::Interpolations.Extrapolation{T,1} = interpolate(cp, mΘ[:,ieps])

        athis = Sthis.itp.knots[1]::Vector{T} ; pop!(athis)                                 # pop last term (there for extrapolation ONLY)
        sthis = Sthis.itp.coefs   ::Vector{T} ; pop!(sthis)

        get_consumption!(cthis, athis, sthis, R, wage, ieps, cp)
        copy!(anext, sthis)

        # loop over next period assets
        MUexp[:] = zero(T)
        for jeps = 1:nEps

            #== Saving policy t+1 ==#
            Snext::Interpolations.Extrapolation{T,1} = interpolate(cp, mΘ′[:,jeps]);

            #== Saving at values ==#
            snext::Vector{T} = Snext[anext]

            #== Consumption next period ==#
            get_consumption!(cnext, anext, snext, R′, wage′, jeps, cp)

            #== Marginal utility next period ==#
            MUnext = dutil(cp, cnext);
            MUexp = MUexp + Π[ieps,jeps] * MUnext;

            #== zero consumption ==#
            cnext[:] = zero(T)
        end

    return inv_dutil(cp, β * R′ * MUexp) , cthis
end
# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

"""
Given a vector construct an interpolation object that does linear interpolation in the asset dimension
"""
function Interpolations.interpolate{T<:Real}(cp::ConsumerProblem, Θ::Vector{T})

    assets  = [ Θ[1]; Θ[1] + cp.asset_knots; 100.0]
    savings = [ 0.0; Θ[2:end]; 0.0 ]

    slope = ( savings[end-1] - savings[end-2] )/( assets[end-1]-assets[end-2] )
    savings[end] = savings[end-1] + ( assets[end]-assets[end-1] ) * slope

    #== Create Interpolant ==#
    S = interpolate( (assets,), savings, Gridded(Linear()) )

    #== extrapolate values outside grid ==#
    S = extrapolate(S, Flat() )
end

function get_consumption!{T<:Real}(C::Vector{T}, asset, savings, R, wage, incomeIndex::Union{AbstractVector{Int64},Integer}, cp::ConsumerProblem)

    unempInd = (incomeIndex .== 1) # BitArray
    effec_wage  = wage * ( cp.μ * unempInd  + (1.0 - cp.τ) * (1 - unempInd) )

    copy!(C, R * asset + effec_wage - savings)
end
function get_consumption(asset, savings, R, wage, incomeIndex::Union{AbstractVector{Int64},Integer}, cp::ConsumerProblem)

    unempInd = (incomeIndex .== 1) # BitArray
    effec_wage  = wage * ( cp.μ * unempInd  + (1.0 - cp.τ) * (1 - unempInd) )

    return R * asset + effec_wage - savings
end

function test_jaco_interp!{T,J<:Real}(out::Vector{T}, Θ::Vector{T}, R::J, wage::J, cp::ConsumerProblem)

    β, Π = cp.β, cp.Π
    nEps, nSavingsPar = cp.nEps, cp.nSavingsPar

    #== Reshape policy coefficients ==#
    mΘ  = reshape(Θ, nSavingsPar, nEps)
    cthis   = Array(T, nSavingsPar)

    # loop over idio shock
    for ieps = 1:nEps
        Sthis::Interpolations.Extrapolation{T,1} = interpolate(cp, mΘ[:,ieps])

        #== CHECK the type ==#
        # println("Type of interpolation object: " , string(typeof(Sthis.itp.coefs[1])) )

        athis = Sthis.itp.knots[1]::Vector{T} ; pop!(athis)                                 # pop last term (there for extrapolation ONLY)
        sthis = Sthis.itp.coefs   ::Vector{T} ; pop!(sthis)

        get_consumption!(cthis, athis, sthis, R, wage, ieps, cp)
        out[(ieps-1)*nSavingsPar+1:ieps*nSavingsPar] = cthis

        cthis[:] = zero(T)
    end

end
