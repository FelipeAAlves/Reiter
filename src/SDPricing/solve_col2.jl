
# **************************************************************************************
#   REVIEW BASICs of CompEcon/BasisMatrices
#
# ======================================================================================
# agrid0 = linspace(0.0.^0.4, 100.0.^0.4, 25).^(1/0.4)
#
# # method two, constructing separately, then calling `Basis` with the two
# y_basis = Basis(Cheb  , 5, -4.0, 4.0)     # Cheb
# a_basis = Basis(Spline, agrid0, 0, 3)
#
# basis = y_basis × a_basis
#
# # Construct state vector (matrix). Note that splidef (called by
# # fundef) adds breakpoints to the original grid we gave it, so let's
# # extract the actual grid points used from the second argument
# S, ( ygrid, agrid) = nodes(basis)
#
# # construct basis matrix and its lu-factorization for very fast inversion
# # NOTE: I am doing this in a round-about way. I could have just done
# #       Φ = BasisMatrix(basis), but doing it this way gives me the direct
# #       representation so I get Φ_y without repeating any calculations
# Φ_tensor = BasisMatrix(basis,Tensor())
#
# #== Compare ways to expand the Tensor ==#
# raw_ind = [collect(1:n) for n in (5, 27)]
# ind = gridmake(raw_ind...)
#
# @time Φ_tensor[1][ind[:,1], :];
# @time kron(ones(27), Φ_tensor[1]);
#
# @time Φ_tensor[2][ind[:,2], :];
# @time kron(Φ_tensor[2], ones(5));
#
# Array( Φ_tensor[2][ind[:,2], :] )
# Array( kron(Φ_tensor[2], ones(5)) )
#
# Φ_direct = BasisMatrix(basis,Direct(), S , 0)
#
# Φ_y  = Array( Φ_direct.vals[1] )
# Φ_y2 = Array( BasisMatrix(y_basis, Expanded()).vals[1] )
# Φ = BasisMatrix(basis,Expanded())
#
# #== FIND coefficients - trick using Expanded or not ==#
# y = repmat(agrid.^2,5) + 2*kron(ones(27), ygrid)
# coeff = funfitxy(basis, S, y)[1];
#
# y1 = Φ.vals[1]*coeff;
# y2 = row_kron(Φ_direct.vals[2], Φ_direct.vals[1])*coeff;
#
# [y y1 y2]
#
# #== Eval at different Points ==#
# agrid_new = linspace(25,75, 10);
# Snew = gridmake(agrid_new, ygrid);
# Φ_new  = BasisMatrix(basis, Expanded(), Snew, 0).vals[1]
# [Φ_new*coeff funeval(coeff, basis, Snew)]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newton_output_file = "C:\\Users\\falves\\Git\\Reiter\\src\\SDPricing\\newton_output.dat"; # []
FF = open(newton_output_file,"w")

"""
  Solve for the firm's value function and policies

"""
function solve_firm_policy!(fcoll::FirmColloc, sol::FirmSolution, w::Float64, Y::Float64=1.0, method::Int64 = 2)

    Φ = fcoll.Φ
    Φ_fac = fcoll.Φ_fac

    # .....................................................................................

    tol = 1e-10
    V    = Φ * sol.coeff
    Vnew = similar(V)

    ii = 1
    err = 1.0
    # global newton_output = open(newton_output_file, "w")
    while err > tol && ii <= 25

        (method==1) && bellman_rhs!(Vnew, fcoll, sol, w, Y)
        (method==2) && bellman_rhs_eval_v!(Vnew, fcoll, sol, w, Y)
        copy!(sol.coeff,Φ_fac\Vnew)

        err = maxabs( vec(V) - vec(Vnew) )
        (mod(ii,100)==0) && @printf("  value function iteration %d, distance %.4f \n", ii, err)

        #== Prepare NEXT iteration ==#
        copy!(V, Vnew)
        fill!(Vnew, 0.0)
        ii += 1
    end
    # close(newton_output)
    return nothing
end

"""
  Computes the RHS of our system of equations for a given guess of collocation
  coefficients. This version creates a function v̂ to use at the evaluations

### Arguments

    - `V`
    - `fcoll`
    - `sol`
"""
function bellman_rhs_eval_v!(V::Array{Float64,2}, fcoll::FirmColloc, sol::FirmSolution, w::Float64, Y::Float64, maximize_pol::Bool = true)

    @getPar __pars

    #== Basis matrices ==#
    # Φ = fcoll.Φ

    #== Indices ==#
    grid_nodes  = fcoll.grid_nodes
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #== Basis matrices ==#
    Φ_z::SparseMatrixCSC{Float64,Int64}   = fcoll.Φ_tensor.vals[2]
    p_basis = fcoll.p_basis

    #== Coefficients for before shock value function ==#
    cv = sol.coeff[:,3]
    function eval_v̂(p̃::Vector{Float64}, z_ind::Vector{Int64}, deriv::Int64 = 0)

        @assert length(p̃)==length(z_ind)
        Φ_p̃ = BasisMatrix( p_basis, Direct(), p̃, deriv).vals[1]
        Φ_eval   = row_kron( Φ_z[ z_ind, :],  Φ_p̃ )
        return Φ_eval * cv
    end

    #== Do th MAX step ==#
    if maximize_pol

        f_val(p̃) = foc_price_adjust_eval_v(p̃, eval_v̂, w, Y)
        f_jac(p̃) = soc_price_adjust_eval_v(p̃, eval_v̂, w, Y)


        ###  INFO:  Newton Iteration     ###
        # println(newton_output, "ITERATIONS\n")
        # pstar = copy(sol.pstar)
        # for it = 1:25
        #     # println(newton_output, f_val(pstar))
        #     @printf(newton_output, " [")
        #     map(x->@printf(newton_output, " % .3e",x), f_val(pstar))
        #     @printf(newton_output, " ]\n")
        #     copy!(pstar, pstar - f_jac(pstar) \ f_val(pstar))
        # end
        # println(newton_output, "------------------------------------\n")

        ###  INFO:  Using Broyden     ###
        (pstar,ind) = broydn(f_val, sol.pstar, [1e-10,1,1], f_jac)

        # (pstar2,ind) = broydn(f_val, sol.pstar, [1e-10,1,1], f_jac)
        # @printf(" DISTANCE: ")
        # map( x -> @printf(" %.3e",x), abs(pstar-pstar2) )
        # @printf("\n")

        # if ind!=0
        #    @printf("  - Problem with foc - residuals: %.3e \n", maxabs(f_val(pstar)) )
        #    pstar = sol.pstar
        # end

        # println(" residuals: \n")
        # map(x -> @printf(" %.3e",x), f_val(pstar))
        # println("\n")

    else
        pstar = log( ϵ/(ϵ-1) * w ./ z_vals )
    end

    #== Compute Value function in case of adjustment ==#
    # value_adjust!(V, pstar, w, Y, cv, p_basis, Φ_z, ind_z_x_z, ind_p_x_z)
    value_adjust_eval_v!(V, pstar,eval_v̂, w, Y)

    #== Compute Value function in case NO adjustment ==#
    value_nadjust_eval_v!(V, eval_v̂, grid_nodes, w, Y)

    #== Adjust/nAdjust decision ==#
    ξstar = ( V[:,1] - V[:,2] ) / w
    ξstar = min(max(ξstar,ξ0), ξbar)

    #== Value beggining period (before taking ξ) ==#
    V[:,3] = ( H(ξstar) .* V[:,1] - w * cond_mean(ξstar) ) + ( 1-H(ξstar) ) .* V[:,2]

    #== UPDATE policies ==#
    copy!(sol.pstar, pstar )
    copy!(sol.ξstar, ξstar )

    return nothing
end

"""
    Get ξ, v at an arbitrary grid for relative prices. Used to back out ξ for the distribution
"""
function get_xi_at{T<:Real}(p_nodes::Vector{Float64}, eval_v̂::Function, pᵃ::Array{Float64}, w::T, Y::T, Y_pri::T=Y, Π_pri::T = 1.0, z_aggr::T=0.0)

    @getPar __pars

    #== Indices ==#
    grid_nodes = gridmake(p_nodes, 1:n_z)             ##  WARN:  index for z_vals     ##

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ##  WARN WARN:  Don't do the maximization          ##
    ##              ONLY check the foc holds at equil  ##
    #== CHECK foc ==#
    # @code_warntype foc_price_adjust_eval_v(sol.pstar, eval_v̂, w, Y, Y_pri, Π_pri, z_aggr)
    # foc_resid = foc_price_adjust_eval_v(sol.pstar, eval_v̂, w, Y, Y_pri, Π_pri, z_aggr)
    # maxabs(foc_resid)<1e-8 || throw(error("foc error of $(maxabs(foc_resid))"))
    # @assert maxabs(foc_resid)<1e-8

    # ##  WARN:  CHECK if code is Dual  ##
    # if eltype(foc_resid)<:ForwardDiff.Dual
    #
    #     ##  Perturb policy  ##
    #     # @code_warntype soc_price_adjust_eval_v(sol.pstar, eval_v̂, w.value, Y.value, Y_pri.value, Π_pri.value)
    #     J = soc_price_adjust_eval_v(sol.pstar, eval_v̂, w.value, Y.value, Y_pri.value, Π_pri.value)
    #     # soc = Float64[ J[i,i].values for i =1:size(J,1) ]
    #
    #     ##  IMPORTANT : Implicit Function Theorem  ##
    #     ##              used just the value for the soc  ##
    #     dpdx  = - foc_resid ./ value(diag(J))
    #     pstar = T[ForwardDiff.Dual(sol.pstar[i], dpdx[i].partials ) for i =1:n_z]
    # else
    #     pstar = copy(sol.pstar)
    # end

    #== Initiliaze V ==#
    V = Array(T, size(grid_nodes,1), 2) ##  NOTE: allowing V to hold type T elements

    #== Compute Value function in case of adjustment ==#
    value_adjust_eval_v!(V, pᵃ, eval_v̂, w, Y, Y_pri, Π_pri, z_aggr)  ##  NOTE: I am using the ENVELOPE here  ##

    #== Compute Value function in case NO adjustment ==#
    value_nadjust_eval_v!(V, eval_v̂, grid_nodes, w, Y, Y_pri, Π_pri, z_aggr)

    #== Adjust/nAdjust decision ==#
    ξstar::Vector{T} = ( V[:,1] - V[:,2] ) / w
    ξstar = min(max(ξstar,ξ0), ξbar)

    #== Value beggining period (before taking ξ) ==#
    Ve::Vector{T} = ( H(ξstar) .* V[:,1] - w * cond_mean(ξstar) ) + ( 1-H(ξstar) ) .* V[:,2]

    ##  WARN WARN:  code     ##

    return ξstar, Ve
end


"""
    Evaluates the foc at vector p̃

### Arguments

    - `p̃`                   vector of rel prices (in log) WARN
    - `v̂_fun::Function`     function with continuation value
    - `w`                   wage
    - `Y`                   output

IMPORTANT: REMEMBER that p̃ stands for log p in the notes
"""
function foc_price_adjust_eval_v{T<:Real}(p̃::Vector{T}, v̂_fun::Function, w::T, Y::T, Y_pri::T=Y, Π_pri::T=1.0, z_aggr::T=0.0)

    @getPar __pars
    # .....................................................................................
    ∂v̂::Vector{T} = v̂_fun( p̃[ind_z_x_z[:,2]] - log(Π_pri), ind_z_x_z[:,1], 1 ) # including effect of inflation/ order of derivative at end
    E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂

    resid = Y * exp(p̃) - ϵ * Y * ( exp(p̃) - w./(exp(z_aggr)*z_vals) ) + β * (Y_pri/Y)^(-σ) * exp( ϵ * p̃) .* ( E∂v̂ )

    return resid
end

function soc_price_adjust_eval_v(p̃::Vector{Float64}, v̂::Function, w::Float64, Y::Float64, Y_pri::Float64=Y, Π_pri::Float64=1.0)

    @getPar __pars

    ∂v̂  = v̂( p̃[ind_z_x_z[:,2]] - log(Π_pri), ind_z_x_z[:,1], 1 )
    ∂²v̂ = v̂( p̃[ind_z_x_z[:,2]] - log(Π_pri), ind_z_x_z[:,1], 2 )

    E∂v̂  = row_kron( eye(n_z) , Π_z ) * ∂v̂
    E∂²v̂ = row_kron( eye(n_z) , Π_z ) * ∂²v̂
    # .....................................................................................

    J = zeros(eltype(∂v̂), n_z, n_z)
    for iz =1:n_z
        J[iz,iz] = Y * (1-ϵ) * exp( p̃[iz]) + β * (Y_pri/Y)^(-σ) *  ϵ * exp( ϵ * p̃[iz]) * E∂v̂[iz]  +
                        β * (Y_pri/Y)^(-σ) * exp( ϵ * p̃[iz]) * E∂²v̂[iz]
    end

    return J

end

function bellman_rhs!(V::Array{Float64,2}, fcoll::FirmColloc, sol::FirmSolution, w::Float64, Y::Float64, max::Bool = true)

    @getPar __pars

    #== Basis matrices ==#
    Φ       = fcoll.Φ
    Φ_z::SparseMatrixCSC{Float64,Int64}   = fcoll.Φ_tensor.vals[2]
    p_basis = fcoll.p_basis

    #== Indices ==#
    grid_nodes  = fcoll.grid_nodes
    ind_z_x_z   = fcoll.ind_z_x_z
    ind_p_x_z   = fcoll.ind_p_x_z

    #== Coefficients for before shock value function ==#
    cv = sol.coeff[:,3]
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if max
        #== FIND optimal p̃ in case of adjustment ==#
        # f!(p̃, rout) = foc_price_adjust!(rout, p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
        # g!(p̃, Jout) = soc_price_adjust!(Jout, p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
        #
        # res = nlsolve(f!, g!, fpol.pstar)
        # if !res.f_converged
        #     @printf("  - Problem with foc - residuals: %.3e \n", res.residual_norm)
        #     pstar = fpol.pstar
        # else
        #     pstar = res.zero
        # end

        f_val(p̃) = foc_price_adjust(p̃, w, Y, cv, p_basis, Φ_z, ind_z_x_z)
        f_jac(p̃) = soc_price_adjust(p̃, w, Y, cv, p_basis, Φ_z, ind_z_x_z)
        # f_val(p̃) = foc_price_adjust_eval_v(p̃, eval_v̂, w, Y, ind_z_x_z)
        # f_jac(p̃) = soc_price_adjust_eval_v(p̃,eval_v̂, w, Y, ind_z_x_z)

        (pstar,ind) = broydn(f_val, sol.pstar, [1e-10,1,1], f_jac) # REVIEW : check options
        if ind!=0
           @printf("  - Problem with foc - residuals: %.3e \n", maxabs(f_val(pstar)) )
           pstar = sol.pstar
        end
    else
        pstar = log( ϵ/(ϵ-1) * w ./ z_vals )
    end

    #== Compute Value function in case of adjustment ==#
    value_adjust!(V, pstar, w, Y, cv, p_basis, Φ_z, ind_z_x_z, ind_p_x_z)
    # value_adjust_eval_v!(V, pstar,eval_v̂, w, Y, ind_z_x_z, ind_p_x_z)

    #== Compute Value function in case NO adjustment ==#
    v̂  = Φ * cv            # continuation value at nodes
    value_nadjust!(V, v̂, grid_nodes, w, Y)

    #== Adjust/nAdjust decision ==#
    ξstar = ( V[:,1] - V[:,2] ) / w
    ξstar[ξstar.>ξbar] = ξbar

    #== Value beggining period (before taking ξ) ==#
    V[:,3] = ( H(ξstar) .* V[:,1] - w * cond_mean(ξstar) ) + ( 1-H(ξstar) ) .* V[:,2]

    #== UPDATE policies ==#
    copy!(sol.pstar, pstar )
    copy!(sol.ξstar, ξstar )

    return nothing
end


function foc_price_adjust!{T<:Real}(resid::Vector{T}, p̃::Vector{T}, w, Y, cv, p_basis, Φ_z, ind_z_x_z)

    @getPar __pars
    # .....................................................................................

    Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
    Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )

    ∂v̂ = Φ * cv
    E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂

    for iz = 1:n_z
        resid[iz] = Y * exp( p̃[iz]) - ϵ * Y * ( exp(p̃[iz]) - w/z_vals[iz]) + β * exp( ϵ * p̃[iz]) * E∂v̂[iz]
    end

    return nothing
end

function foc_price_adjust{T<:Real}(p̃::Vector{T}, w, Y, cv, p_basis, Φ_z, ind_z_x_z)

    @getPar __pars
    # .....................................................................................

    Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
    Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )

    ∂v̂ = Φ * cv
    E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂

    resid = Y* exp(p̃) - ϵ * Y * ( exp(p̃) - w./z_vals) + β * exp(ϵ * p̃) .* E∂v̂

    return resid
end


function soc_price_adjust!{T<:Real}(J::Matrix{T}, p̃::Vector{T}, w, Y, cv, p_basis, Φ_z, ind_z_x_z)

    @getPar __pars

    Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
    Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
    ∂v̂ = Φ * cv
    # .....................................................................................
    Φ_p̃_deriv2 = BasisMatrix( p_basis, Direct(), p̃, 2).vals[1]
    Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
    ∂²v̂ = Φ * cv

    E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂
    E∂²v̂ = row_kron( eye(n_z) , Π_z ) * ∂²v̂
    # .....................................................................................

    fill!(J,zero(T))
    for iz =1:n_z
        J[iz,iz] = Y * (1-ϵ) * exp( p̃[iz]) + β * ϵ * exp( ϵ * p̃[iz]) * E∂v̂[iz] +
                        β * exp( ϵ * p̃[iz]) * E∂²v̂[iz]
    end

    return nothing
end

function soc_price_adjust{T<:Real}(p̃::Vector{T}, w, Y, cv, p_basis, Φ_z, ind_z_x_z)

    @getPar __pars

    Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
    Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
    ∂v̂ = Φ * cv
    # .....................................................................................
    Φ_p̃_deriv2 = BasisMatrix( p_basis, Direct(), p̃, 2).vals[1]
    Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
    ∂²v̂ = Φ * cv

    E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂
    E∂²v̂ = row_kron( eye(n_z) , Π_z ) * ∂²v̂
    # .....................................................................................

    J = zeros(T, n_z, n_z)
    for iz =1:n_z
        J[iz,iz] = Y * (1-ϵ) * exp( p̃[iz]) + β * ϵ * exp( ϵ * p̃[iz]) * E∂v̂[iz] +
                        β * exp( ϵ * p̃[iz]) * E∂²v̂[iz]
    end

    return J
end



# function foc_price_adjust(p̃, z_vals, w,  cv, fcoll, p_basis, Φ_z, ind_z_x_z)
#
#     n_z = length(z_vals)
#     # .....................................................................................
#
#     Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
#     Φ = row_kron( Φ_p̃_deriv[ ind_z_x_z[:,2], :], Φ_z[ ind_z_x_z[:,1], :] )
#     ∂v̂ = Φ * cv
#
#     E∂v̂ = row_kron( eye(n_z) , fcoll.Π_z ) * ∂v̂
#
#     return exp(-ϵ * p̃) * Y - ϵ * ( exp(p̃) - w/z_vals ) * exp( -(ϵ + 1.0) * p̃) * Y + β * E∂v̂
# end


function value_adjust!(V::Array{Float64,2}, p̃::Vector{Float64}, w, Y, cv, p_basis, Φ_z, ind_z_x_z, ind_p_x_z)

    @getPar __pars
    n_p̃ = length(p̃)
    # .....................................................................................

    #== Construct Basis at p̃ ==#
    Φ_p̃ = BasisMatrix( p_basis, Direct(), p̃, 0).vals[1]
    Φ   = row_kron( Φ_z[ ind_z_x_z[:,1], :],  Φ_p̃[ ind_z_x_z[:,2], :] )

    #== Find the Expected continuation value ==#
    v̂  = Φ * cv
    Ev̂ = row_kron( eye(n_z), Π_z ) * v̂

    val_adjust = profit(p̃, z_vals, w, Y) + β * Ev̂

    V[:,1] = val_adjust[ ind_p_x_z[:,2] ]

    return nothing
end


function value_adjust!{T<:Real}(V::Array{T,2}, p̃::Vector{Float64}, w::T, Y::T, cv, p_basis, Φ_z, ind_p_x_z)

    @getPar __pars
    n_p̃ = length(p̃)
    # .....................................................................................

    Φ_p̃ = BasisMatrix( p_basis, Direct(), p̃, 0).vals[1]
    Ev̂  = zeros(n_p̃)
    zind = zeros(n_p̃)
    for iz = 1:n_z
        fill!(zind,iz)
        Ev̂ = Ev̂ + Π_z[:,iz] .* ( row_kron( Φ_z[ zind, :] , Φ_p̃ ) * cv)
    end

    val_adjust = profit(p̃, z_vals, w, Y) + β * Ev̂

    V[:,1] = val_adjust[ind_p_x_z[:,2]]
end

"""
p̃ in log units
"""
function value_adjust_eval_v!{J<:Real,T<:Real}(V::Array{T,2}, p̃::Vector{J}, v̂::Function, w::T, Y::T, Y_pri::T=Y, Π_pri::T=1.0, z_aggr::T=0.0)

    @getPar __pars
    n_p = div(size(V,1), n_z) ##  NOTE:  size of price index  ##
    # .....................................................................................

    #== Find the Expected continuation value ==#
    Ev̂::Vector{T} = row_kron( eye(n_z), Π_z ) * v̂( p̃[ind_z_x_z[:,2]] - log(Π_pri), ind_z_x_z[:,1] )

    val_adjust = profit(p̃, z_vals, w, Y, z_aggr) + β * (Y_pri/Y)^(-σ) * Ev̂

    V[:,1] = val_adjust[ kron(1:n_z,ones(Int64,n_p)) ]

    return nothing
end

function value_nadjust!(V::Array{Float64}, v̂::Array{Float64}, grid, w, Y)

    @getPar __pars

    #== Find the continuation expected value ==#
    Ev̂ = v2Ev(v̂, Π_z)

    p̃m1     = grid[:,1]
    z_ind::Vector{Int64} = grid[:,2]

    V[:,2] = profit(p̃m1, z_vals[z_ind], w, Y) + β * Ev̂

    return nothing
end

function value_nadjust_eval_v!{T<:Real}(V::Array{T}, v̂::Function, grid::Array{Float64}, w::T, Y::T, Y_pri::T=Y, Π_pri::T=1.0, z_aggr::T=0.0)

    @getPar __pars

    xm1 = grid[:,1]
    z_ind::Vector{Int64} = grid[:,2]

    #== Find the continuation expected value ==#
    Ev̂::Vector{T} = v2Ev( v̂(xm1 - log(Π_pri), z_ind), Π_z )

    V[:,2] = profit(xm1, z_vals[z_ind], w, Y, z_aggr) + β * (Y_pri/Y)^(-σ) * Ev̂

    return nothing
end

function profit(p̃, z_nodes, w, Y, z_aggr)
    @getPar __pars

    return ( exp(p̃) - w./ (exp(z_aggr)*z_nodes) ) .* ( exp( -ϵ * p̃)  ) * Y
end

function v2Ev(v,Π_z)

    n_z = size(Π_z,1)
    n_p = div(length(v),n_z )

    return vec( reshape(v, n_p, n_z) * Π_z' )

end

function v2Ev!(Ev, v, Π_z)

    n_z = size(Π_z,1)
    n_p = div(length(v),n_z )

    copy!(Ev, vec(reshape(v, n_p, n_z) * Π_z') )

end
# function profit_deriv(p̃, z_vals,w, Y, fcoll)
#
#     return ( exp(p̃) - w./ z_vals ) .* ( exp(p̃) .^( -ϵ ) ) * Y
# end
# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------
