"""
Computes the residual of asset market using histogram approximation of distribution.
The law of motion for histogram is done as in Young (2010);

### INPUT

"""
function stst_histogram_resid(w::Float64, Y::Float64, ss_histogram::StstHistogram, fcoll::FirmColloc, sol::FirmSolution, update_ss::Bool=false)

    @getPar __pars
    hist_nodes, (p_hist_nodes, z_hist_nodes) = nodes(ss_histogram)

    # .....................................................................................

    #=======================================================#
    ### Compute Policy function for set of prices         ###
    #=======================================================#
    @printf("   - Inner loop (Policy) \n")

    #== Save and reshape solution ==#
    @time solve_firm_policy!(fcoll, sol, w, Y)

    #== TO compute eval_v̂ ==#
    Φ_z::SparseMatrixCSC{Float64,Int64}   = fcoll.Φ_tensor.vals[2]
    p_basis = fcoll.p_basis
    cv = sol.coeff[:,3]
    function eval_v̂(x::Vector{Float64}, z_ind::Vector{Int64}, deriv::Int64 = 0)

        @assert length(x)==length(z_ind)
        Φ_p̃ = BasisMatrix( p_basis, Direct(), x, deriv).vals[1]
        Φ_eval   = row_kron( Φ_z[ z_ind, :],  Φ_p̃ )
        return Φ_eval * cv
    end
    ξstar_distr, _ = get_xi_at(p_hist_nodes, eval_v̂, sol.pstar, w, Y)

    #=================================================================#
    ### Compute Stationary Distribution from decision rules         ###
    #=================================================================#

    #== Compute Transpose Transition matrix ==#
    Π = endperiod_transition( p_hist_nodes, sol.pstar, ξstar_distr)

    #== Invariant distribution ==#
    # @time _, eigenv = eigs(Π, nev = 1,which=:LM);
    # eigenv = real( squeeze(eigenv,2) )
    # vHistogram = eigenv ./ sum(eigenv);

    vHistogram   = zeros(size(Π,1));
    vHistogram[1] = 1.0;
    vHistogram_t1= similar(vHistogram);
    i = 1
    err_hist = 1.0
    while err_hist>1e-10 && i<10000
        copy!(vHistogram_t1,Π*vHistogram)
        err_hist = maxabs(vHistogram - vHistogram_t1)
        copy!(vHistogram,vHistogram_t1)
        i += 1
    end

    #== CHECK all elements are positive ==#
    minimum(vHistogram) > - 1e-10 || throw(error("Negative histogram weight"))
    vHistogram[vHistogram .< 1e-10] = 0.0
    vHistogram = vHistogram ./ sum(vHistogram)

    #== CHECK consistency ==#
    err_invDist = maximum( abs( Π*vHistogram - vHistogram) );
    @printf("   - Invariant distribution Converged to %.2e \n", err_invDist)
    err_invDist<1e-9 || throw(error("Invariante distribution is incorrect"))

    # vEndHistogram = Π_in * vHistogram
    # mEndHistogram = reshape(vEndHistogram, length(p_hist_nodes), length(z_hist_nodes))

    #===================================================#
    ###   COMPUTE residual on equilibrium condition   ###
    #===================================================#
    ##  WARN WARN:  END-of period prices = beginning of period price due     ##
    ##              to NO inflation at SS                                    ##
    Π_in = inperiod_transition(p_hist_nodes, sol.pstar, ξstar_distr ; update_idio=false)
    vHistogram_end = Π_in*vHistogram

    #== pricing function ==#
    n_hist_p = div(length(vHistogram),n_z)
    p_histogram_end = sum(reshape(vHistogram_end, n_hist_p, n_z),2)
    p_histogram_end = squeeze(p_histogram_end,2)
    pricing_fnc_resid = pricing_fnc(hist_nodes, sol.pstar, ξstar_distr, vHistogram, p_histogram_end)[1]

    #== Labor market ==#
    labor_mkt_resid   = resid_labor(hist_nodes, ξstar_distr, vHistogram, vHistogram_end, Y)

    @printf("   ∫p^(1-ϵ)       : %8.6f\n", pricing_fnc_resid)
    @printf("   N^* - x        : %8.6f\n", labor_mkt_resid)
    println("---------------------------------------")
    #
    if update_ss

        #== Update ss_histogram variable ==#
        copy!(ss_histogram.vHistogram, vHistogram)
        copy!(ss_histogram.mHistogram, vHistogram)

        ##  WARN WARN:  expect grid_nodes     ##
        x     = fcoll.grid_nodes[:,1]
        z_ind::Vector{Int64} = fcoll.grid_nodes[:,2]
        Ve_bellman = eval_v̂(x,z_ind)

        ss_histogram.xstst = [vHistogram[2:end]; 0.0; 0.0]                      # [vHistogram, Z, ϵ]
        ss_histogram.ystst = [Y;Nstst; 1/β-1.0; 1.0; w; Ve_bellman; sol.pstar]  # [Y_t, N_t, i_t, Π_t, w_t, Ve, pᵃ]
        ss_histogram.χ = Y^(-σ)*w/(Nstst^(1/ϕ)); # set the χ
        return Void
    else
        return [pricing_fnc_resid; labor_mkt_resid]
    end
end


"""
    Price fnc definition

"""
function pricing_fnc{T<:Real}(hist_nodes::Array{Float64}, pᵃ::Vector{T}, ξstar::Vector{T}, vHistogram_begin::Vector{T}, p_histogram_end::Vector{T}, verify::Bool = true)

    @getPar __pars
    p_nodes = hist_nodes[:,1]
    n_hist_p = div(length(p_nodes),n_z)

    #== end of period distribution over price ==#
    pricing_fnc_resid1::T = 1.0 - dot( exp( (1-ϵ)*p_nodes[1:n_hist_p]) , p_histogram_end )
    #== beginning of period distribution ==#
    pricing_fnc_resid2::T = 1.0 - dot( H(ξstar) .* kron(exp( (1-ϵ)*pᵃ), ones(n_hist_p)) +  (1-H(ξstar)).* exp( (1-ϵ)*p_nodes ), vHistogram_begin)
    if verify
        if abs(pricing_fnc_resid1 - pricing_fnc_resid2)>1e-8
            @printf("Different results in pricing equilibrium condition\n")
            @printf(" 01: %8.6f\n", pricing_fnc_resid1 )
            @printf(" 02: %8.6f\n", pricing_fnc_resid2 )
        end
        return pricing_fnc_resid1,pricing_fnc_resid2
    else
        return pricing_fnc_resid1
    end

end

"""
Residual on the labor market

    N = labor_demand + cost_of_changing_prices

    - First  term depends on end_period distribution of prices/technology
    - second term depends on begin_period distribution of prices/technology
"""
function resid_labor{T<:Real}(hist_nodes::Array{Float64}, ξstar::Vector{T}, vHistogram_begin::Vector{T}, vHistogram_end::Vector{T}, Y::T, z_aggr::T=0.0)
    @getPar __pars

    p_nodes, z_ind_nodes::Vector{Int64} = hist_nodes[:,1], hist_nodes[:,2]
    resid::T = Nstst - dot( Y * exp(-ϵ * p_nodes ) ./(exp(z_aggr) * z_vals[z_ind_nodes]) , vHistogram_end) - dot( cond_mean(ξstar) , vHistogram_begin)

    return resid
end

"""
Compute the TRANSPOSE of transition matrix BETWEEN-periods from policies

IMPORTANT
    - distribution starts over
        * t-1 END period relative prices p̃₁
        * t period idio shocks  a_t
    - evolves one period. That involves 3 adjustments
        * scale of inflation at t           - x = p_{t-1}/P_t
        * choices from beginning period t   - p̃(p_{t-1}/P_t)
        * evol of idio shocks               - a_t --> a_{t+1} (OPTION update_idio)

These steps are separated into tow different functions

    * Πadj_transition     :  inflation adjustment
    * inperiod_transition :  last two adjustments
"""
function endperiod_transition{T<:Real}( p_nodes::Array{Float64}, pstar::Vector{T}, ξstar::Vector{T}, Π::T=1.0; update_idio::Bool=true )


    return inperiod_transition(p_nodes, pstar, ξstar; update_idio=update_idio) * Πadj_transition(p_nodes, Π)
    # return Ir, Ic, Val
end
"""
Compute the TRANSPOSE of transition matrix IN-the-period from policies

"""
function inperiod_transition{T<:Real}(p_nodes::Array{Float64}, pstar::Vector{T}, ξstar::Vector{T} ; update_idio::Bool=true)

    @getPar __pars

    n_p = length(p_nodes)
    I = eye(n_z)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ###   NOTE:  real states today due to inflation   ###
    p_nodes = p_nodes

    #== Initialize useful vector ==#
    ic1::Vector{Int64} = collect( 1:length(p_nodes) )
    ic2::Vector{Int64} = kron(collect( 1:length(p_nodes) ), ones(Int64,2))

    #== indices ==#
    i_inaction_to = ic1
    i_action_to   = Array(Int64, 2, n_p)

    #== weights ==#
    vv_action   = Array(T, 2, n_p)

    #== Initialize Vectors for sparse matrix ==#
    Ic  = Array(Int64,0)
    Ir  = Array(Int64,0)
    Val = Array(T,0)         # should contain Reals

    ###  NOTE:  loop over shocks today state-TODAY     ###
    for iz = 1:n_z

        #== states today start at ==#
        ioff   = (iz-1) * n_p
        ind = ioff+1:iz*n_p

        #== p_inac at each relevant node ==#
        ξ = repmat( H( ξstar[ind] ), 1, 2)

        #== Inaction case ==#
        # linear_trans!(i_inaction_to, vv_inaction, p_nodes, p_nodes )

        #== Action case ==#
        linear_trans!(i_action_to, vv_action, p_nodes, ones(n_p) * pstar[iz] )

        for jz = 1:n_z ###  NOTE: shock tomorrow  ###

            p_trans_exo  = (update_idio) * Π_z[iz, jz] + (~update_idio) * I[iz,jz]
            joff   = (jz-1) * n_p

            if p_trans_exo>0.0
                #== Inaction case ==#
                append!(Ic, ic1 + ioff)
                append!(Ir, i_inaction_to + joff)
                append!(Val, (1.0 - ξ[:,1]) * p_trans_exo )

                #== Action case ==#
                append!(Ic, ic2 + ioff)
                append!(Ir, i_action_to + joff)
                append!(Val, ( transpose(ξ) .* vv_action) * p_trans_exo )
            end
        end
    end
    nΠ = n_p * n_z

    return sparse(Ir, Ic, Val, nΠ, nΠ)
end

function Πadj_transition{T<:Real}(p_nodes::Array{Float64}, Π::T)

    @getPar __pars
    n_p = length(p_nodes)

    #== Initialize useful vector ==#
    ic::Vector{Int64} = kron(collect( 1:length(p_nodes) ), ones(Int64,2))
    ito = Array(Int64, 2, n_p)
    vv = Array(T, 2, n_p)

    #== Initialize Vectors for sparse matrix ==#
    Ic  = Array(Int64,0)
    Ir  = Array(Int64,0)
    Val = Array(T,0)         # should contain Reals

    ##  NOTE:  loop over shocks today     ##
    for iz = 1:n_z

        #== states today start at ==#
        ioff   = (iz-1) * n_p

        #== transition due to inflation ==#
        linear_trans!(ito, vv, p_nodes, p_nodes/Π)

        ###  WARN:  no idio tranistion on this step     ###

        #== Append on indices ==#
        append!(Ic, ic + ioff)
        append!(Ir, ito + ioff)
        append!(Val, vv )
    end
    nΠ = n_p * n_z

    return sparse(Ir, Ic, Val, nΠ, nΠ)
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
    nn = length(grid)

    for (i,value) in enumerate(pol)
        # println("($i,$value)")
        ##  NOTE: 1st value ≥ in grid
        ipos   = min( max( searchsortedfirst(grid, value), 2), nn)

        iTo[i]   = ipos
        pHi = ( value - grid[ipos-1] ) / ( grid[ipos] - grid[ipos-1] )
        pHigh[i] = min(max(pHi,0.0), 1.0)

    end

    #== Complete indices ==#
    iFr  = [iFr iFr ]
    iTo  = [iTo (iTo-1) ]
    prob = [pHigh (1.0 - pHigh)]

    return iFr, iTo, prob
end
"""
Find the index going to and the weights on i, i+1
"""
function linear_trans!{T<:Real}(iTo::Array{Int64}, prob::Array{T}, grid::Vector{Float64}, pol::Vector{T})

    # @assert size(iTo,1)==length(pol)

    if size(iTo,2)==1 # return only p_high
        for (i,value) in enumerate(pol)
            # println("($i,$value)")

            ##  WARN:  1st value ≥ in grid     ##
            ipos   = min( max( searchsortedfirst(grid, value), 2), length(grid))

            iTo[i]    = ipos

            ## WARN : Adjust for values outside Grid  ##
            pp::T = ( value - grid[ipos-1] ) / ( grid[ipos] - grid[ipos-1] )
            pp    = min( max(pp,0.0), 1.0 )
            prob[i] = pp
        end
    elseif size(iTo,2)>1 # return [p_high 1.0-p_high]
        for (i,value) in enumerate(pol)
            # println("($i,$value)")

            ##  WARN:  1st value ≥ in grid     ##
            ipos   = min( max( searchsortedfirst(grid, value), 2), length(grid))
            iTo[:,i]    = [ipos;ipos-1]     # two row vectors

            ## WARN : Adjust for values outside Grid  ##
            pp::T = ( value - grid[ipos-1] ) / ( grid[ipos] - grid[ipos-1] )
            pp    = min( max(pp,0.0), 1.0 )

            prob[:,i] = [pp; 1.0-pp]
        end
    end

    return nothing

end


###  Interpolations test  ###
# **************************************************************************************
#== Policies at state today (scaled by inflation) ==#
# ξ_itp::Interpolations.GriddedInterpolation{T,1,T,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1}},0} =
#         interpolate( (p_nodes,), ξstar[ind], Gridded(Linear()) )
# ξ_itp::Interpolations.Extrapolation{T,1,Interpolations.GriddedInterpolation{T,1,T,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1}},0},Interpolations.Gridded{Interpolations.Linear},Interpolations.OnGrid,Interpolations.Flat} =
#         extrapolate( ξ_itp, Flat() )
###  INFO : test for this below  ###

# pnodes = [0.5;1.0;1.5]
# ξstar  = [ForwardDiff.Dual(1.0,(0.0,1.0,0.0,0.0));
#           ForwardDiff.Dual(1.5,(0.0,0.0,1.0,0.0));
#           ForwardDiff.Dual(2.0,(0.0,0.0,0.0,1.0))]
# ξstar  = [1.0;1.5;2.0];
#
# ξ_itp = interpolate( (pnodes,), ξstar, Gridded(Linear()) )
# ξ_itp[0.00]
# ξ_itp[0.75]
# ξ_itp[1.2]
# ξ_itp = extrapolate( ξ_itp, Flat() )
# ξ_itp[0.00]
# ξ_itp[0.75]
# ξ_itp[1.2]
#
# p_Πadj_nodes = ForwardDiff.Dual(0.75,(-1.0,0.0,0.0,0.0))
# ξ_itp[0.75]
# ξ_itp[p_Πadj_nodes]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
