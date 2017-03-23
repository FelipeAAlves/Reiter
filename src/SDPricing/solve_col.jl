#
# # **************************************************************************************
# #   REVIEW BASICs of CompEcon/BasisMatrices
# #
# # ======================================================================================
# # agrid0 = linspace(0.0.^0.4, 100.0.^0.4, 25).^(1/0.4)
# #
# # # method two, constructing separately, then calling `Basis` with the two
# # y_basis = Basis(Cheb  , 5, -4.0, 4.0)     # Cheb
# # a_basis = Basis(Spline, agrid0, 0, 3)
# #
# # basis = y_basis × a_basis
# #
# # # Construct state vector (matrix). Note that splidef (called by
# # # fundef) adds breakpoints to the original grid we gave it, so let's
# # # extract the actual grid points used from the second argument
# # S, ( ygrid, agrid) = nodes(basis)
# #
# # # construct basis matrix and its lu-factorization for very fast inversion
# # # NOTE: I am doing this in a round-about way. I could have just done
# # #       Φ = BasisMatrix(basis), but doing it this way gives me the direct
# # #       representation so I get Φ_y without repeating any calculations
# # Φ_tensor = BasisMatrix(basis,Tensor())
# #
# # #== Compare ways to expand the Tensor ==#
# # raw_ind = [collect(1:n) for n in (5, 27)]
# # ind = gridmake(raw_ind...)
# #
# # @time Φ_tensor[1][ind[:,1], :];
# # @time kron(ones(27), Φ_tensor[1]);
# #
# # @time Φ_tensor[2][ind[:,2], :];
# # @time kron(Φ_tensor[2], ones(5));
# #
# # Array( Φ_tensor[2][ind[:,2], :] )
# # Array( kron(Φ_tensor[2], ones(5)) )
# #
# # Φ_direct = BasisMatrix(basis,Direct(), S , 0)
# #
# # Φ_y  = Array( Φ_direct.vals[1] )
# # Φ_y2 = Array( BasisMatrix(y_basis, Expanded()).vals[1] )
# # Φ = BasisMatrix(basis,Expanded())
# #
# # #== FIND coefficients - trick using Expanded or not ==#
# # y = repmat(agrid.^2,5) + 2*kron(ones(27), ygrid)
# # coeff = funfitxy(basis, S, y)[1];
# #
# # y1 = Φ.vals[1]*coeff;
# # y2 = row_kron(Φ_direct.vals[2], Φ_direct.vals[1])*coeff;
# #
# # [y y1 y2]
# #
# # #== Eval at different Points ==#
# # agrid_new = linspace(25,75, 10);
# # Snew = gridmake(agrid_new, ygrid);
# # Φ_new  = BasisMatrix(basis, Expanded(), Snew, 0).vals[1]
# # [Φ_new*coeff funeval(coeff, basis, Snew)]
#
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# """
#     Solve for the firm's value function and policies
# """
# function solve_firm_policy!(w::Float64, fcoll::FirmColloc, sol::FirmSolution)
#
#     Φ = fcoll.Φ
#     Φ_fac = fcoll.Φ_fac
#     # .....................................................................................
#
#     tol = 1e-6
#     V    = Φ * sol.coeff
#     Vnew = similar(V)
#
#     ii = 1
#     err = 1.0
#     while err > tol && ii <= 1000
#
#         bellman_rhs!(Vnew, w, fcoll, sol)
#         copy!(sol.coeff,Φ_fac\Vnew)
#
#         err = maxabs( vec(V) - vec(Vnew) )
#         (mod(ii,10)==0) && @printf("  value function iteration %d, distance %.4f \n", ii, err)
#
#         #== Prepare NEXT iteration ==#
#         copy!(V, Vnew)
#         fill!(Vnew, 0.0)
#         ii += 1
#     end
#
#     return Void
# end
#
# """
# Computes the RHS of our system of equations for a given guess of collocation
# coefficients.
#
# ## Arguments
#
# - `V`
#
# """
# function bellman_rhs!(V::Array{Float64,2}, w::Float64, fcoll::FirmColloc, sol::FirmSolution, max::Bool = true)
#
#     @getPar __pars
#
#     #== Basis matrices ==#
#     Φ       = fcoll.Φ
#     Φ_z   = fcoll.Φ_tensor.vals[2]
#     p_basis = fcoll.p_basis
#
#     #== Indices ==#
#     grid_nodes = fcoll.grid_nodes
#     ind_z_x_z = fcoll.ind_z_x_z
#     ind_p_x_z = fcoll.ind_p_x_z
#
#     #== Coefficients for before shock value function ==#
#     cv = sol.coeff[:,3]
#     v̂ = Φ * cv     # continuation value
#     # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#     if max
#         #== FIND optimal p̃ in case of adjustment ==#
#         # f!(p̃, rout) = foc_price_adjust!(rout, p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
#         # g!(p̃, Jout) = soc_price_adjust!(Jout, p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
#         #
#         # res = nlsolve(f!, g!, fpol.pstar)
#         # if !res.f_converged
#         #     @printf("  - Problem with foc - residuals: %.3e \n", res.residual_norm)
#         #     pstar = fpol.pstar
#         # else
#         #     pstar = res.zero
#         # end
#
#         f_val(p̃) = foc_price_adjust(p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
#         f_jac(p̃) = soc_price_adjust(p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
#
#         (pstar,ind) = broydn(f_val, sol.pstar, [1e-10,1,1], f_jac)
#         if ind!=0
#            @printf("  - Problem with foc - residuals: %.3e \n", maxabs(f_val(pstar)) )
#            pstar = sol.pstar
#         end
#     else
#         pstar = log( ϵ/(ϵ-1) * w ./ z_vals )
#     end
#
#     #== Compute Value function in case of adjustment ==#
#     value_adjust!(V, pstar, w, cv, p_basis, Φ_z, ind_z_x_z, ind_p_x_z)
#
#     #== Compute Value function in case NO adjustment ==#
#     value_nadjust!(V, grid_nodes, w, v̂)
#
#     #== Adjust/nAdjust decision ==#
#     ξstar = ( V[:,1] - V[:,2] ) / w
#     ξstar[ξstar.>ξbar] = ξbar
#
#     #== Value beggining period (before taking ξ) ==#
#     V[:,3] = ( H(ξstar) .* V[:,1] - w * cond_mean(ξstar) ) + ( 1-H(ξstar) ) .* V[:,2]
#
#     #== UPDATE policies ==#
#     copy!(sol.pstar, pstar )
#     copy!(sol.ξstar, ξstar )
#
#     return Void
# end
#
# """
# Get policies at an arbitrary grid
# """
# # REVIEW : Check if fnc works as intended
# function get_xi_at(p_nodes::Vector{Float64}, z_nodes::Vector{Float64}, w::Float64, fcoll::FirmColloc, sol::FirmSolution)
#
#     @getPar __pars
#
#     #== Coefficients for continuation value ==#
#     cv = sol.coeff[:,3]
#
#     #== Basis matrices ==#
#     basis   = fcoll.basis
#     Φ_z   = fcoll.Φ_tensor.vals[2]
#     p_basis = fcoll.p_basis
#
#     #== Indices ==#
#     ind_z_x_z  = fcoll.ind_z_x_z
#     ind_p_x_z  = gridmake(1:length(p_nodes), 1:length(z_nodes))  ##  NOTE: dif indice here  ##
#     grid_nodes = gridmake(p_nodes,z_nodes)
#     # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#     f_val(p̃) = foc_price_adjust(p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
#     f_jac(p̃) = soc_price_adjust(p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
#
#     (pstar,ind) = broydn(f_val, sol.pstar, [1e-10,1,1], f_jac)
#     if ind!=0
#         @printf("  - Problem with foc - residuals: %.3e \n", maxabs(f_val(pzero)) )
#         return
#     end
#
#     #== Initiliaze V ==#
#     V = zeros( size(grid_nodes,1), 2 )
#
#     #== Compute Value function in case of adjustment ==#
#     value_adjust!(V, pstar, w, cv, p_basis, Φ_z, ind_z_x_z, ind_p_x_z)
#
#     #== Compute Value function in case NO adjustment ==#
#     Φ_mat = BasisMatrix(basis, Expanded(), grid_nodes).vals[1]
#     v̂     = Φ_mat * cv
#     value_nadjust!(V, grid_nodes, w, v̂)
#
#     #== Adjust/nAdjust decision ==#
#     ξstar = ( V[:,1] - V[:,2] ) / w
#     ξstar[ξstar.>ξbar] = ξbar
#
#     #== UPDATE policies ==#
#     sol.pol_at[1] = (p_nodes, pstar, ξstar)
#
#     return nothing
# end
#
#
# function foc_price_adjust!{T<:Real}(resid::Vector{T}, p̃::Vector{T}, w, cv, p_basis, Φ_z, ind_z_x_z)
#
#     @getPar __pars
#     # .....................................................................................
#
#     Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
#     Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
#
#     ∂v̂ = Φ * cv
#     E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂
#
#     for iz = 1:n_z
#         resid[iz] = exp( p̃[iz]) - ϵ * ( exp(p̃[iz]) - w/z_vals[iz]) + β * exp( ϵ * p̃[iz]) * E∂v̂[iz]
#     end
#
#     return Void
# end
#
# function foc_price_adjust{T<:Real}(p̃::Vector{T}, w, cv, p_basis, Φ_z, ind_z_x_z)
#
#     @getPar __pars
#     # .....................................................................................
#
#     Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
#     Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
#
#     ∂v̂ = Φ * cv
#     E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂
#
#     resid = exp( p̃) - ϵ * ( exp(p̃) - w./z_vals) + β * exp(ϵ * p̃) .* E∂v̂
#
#     return resid
# end
#
# function soc_price_adjust!{T<:Real}(J::Matrix{T}, p̃::Vector{T}, w, cv, p_basis, Φ_z, ind_z_x_z)
#
#     @getPar __pars
#
#     Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
#     Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
#     ∂v̂ = Φ * cv
#     # .....................................................................................
#     Φ_p̃_deriv2 = BasisMatrix( p_basis, Direct(), p̃, 2).vals[1]
#     Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
#     ∂²v̂ = Φ * cv
#
#     E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂
#     E∂²v̂ = row_kron( eye(n_z) , Π_z ) * ∂²v̂
#     # .....................................................................................
#
#     fill!(J,zero(T))
#     for iz =1:n_z
#         J[iz,iz] = (1-ϵ) * exp( p̃[iz]) + β * ϵ * exp( ϵ * p̃[iz]) * E∂v̂[iz] +
#                         β * exp( ϵ * p̃[iz]) * E∂²v̂[iz]
#     end
#
#     return Void
# end
#
# function soc_price_adjust{T<:Real}(p̃::Vector{T}, w, cv, p_basis, Φ_z, ind_z_x_z)
#
#     @getPar __pars
#
#     Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
#     Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
#     ∂v̂ = Φ * cv
#     # .....................................................................................
#     Φ_p̃_deriv2 = BasisMatrix( p_basis, Direct(), p̃, 2).vals[1]
#     Φ = row_kron( Φ_z[ ind_z_x_z[:,1], :], Φ_p̃_deriv[ ind_z_x_z[:,2], :] )
#     ∂²v̂ = Φ * cv
#
#     E∂v̂ = row_kron( eye(n_z) , Π_z ) * ∂v̂
#     E∂²v̂ = row_kron( eye(n_z) , Π_z ) * ∂²v̂
#     # .....................................................................................
#
#     J = zeros(T, n_z, n_z)
#     for iz =1:n_z
#         J[iz,iz] = (1-ϵ) * exp( p̃[iz]) + β * ϵ * exp( ϵ * p̃[iz]) * E∂v̂[iz] +
#                         β * exp( ϵ * p̃[iz]) * E∂²v̂[iz]
#     end
#
#     return J
# end
#
#
# # function foc_price_adjust(p̃, z_vals, w,  cv, fcoll, p_basis, Φ_z, ind_z_x_z)
# #
# #     n_z = length(z_vals)
# #     # .....................................................................................
# #
# #     Φ_p̃_deriv = BasisMatrix( p_basis, Direct(), p̃, 1).vals[1]
# #     Φ = row_kron( Φ_p̃_deriv[ ind_z_x_z[:,2], :], Φ_z[ ind_z_x_z[:,1], :] )
# #     ∂v̂ = Φ * cv
# #
# #     E∂v̂ = row_kron( eye(n_z) , fcoll.Π_z ) * ∂v̂
# #
# #     return exp(-ϵ * p̃) * Y - ϵ * ( exp(p̃) - w/z_vals ) * exp( -(ϵ + 1.0) * p̃) * Y + β * E∂v̂
# # end
#
# function value_adjust!(V::Array{Float64,2}, p̃::Vector{Float64}, w, cv, p_basis, Φ_z, ind_z_x_z, ind_p_x_z)
#
#     @getPar __pars
#     n_p̃ = length(p̃)
#     # .....................................................................................
#
#     #== Construct Basis at p̃ ==#
#     Φ_p̃ = BasisMatrix( p_basis, Direct(), p̃, 0).vals[1]
#     Φ   = row_kron( Φ_z[ ind_z_x_z[:,1], :],  Φ_p̃[ ind_z_x_z[:,2], :] )
#
#     #== Find the Expected continuation value ==#
#     v̂  = Φ * cv
#     Ev̂ = row_kron( eye(n_z), Π_z ) * v̂
#
#     val_adjust = profit(p̃, z_vals, w) + β * Ev̂
#
#     V[:,1] = val_adjust[ind_p_x_z[:,2]]
# end
#
# function value_adjust!(V::Array{Float64,2}, p̃::Vector{Float64}, w, cv, p_basis, Φ_z, ind_p_x_z)
#
#     @getPar __pars
#     n_p̃ = length(p̃)
#     # .....................................................................................
#
#     Φ_p̃ = BasisMatrix( p_basis, Direct(), p̃, 0).vals[1]
#     Ev̂  = zeros(n_p̃)
#     zind = zeros(n_p̃)
#     for iz = 1:n_z
#         fill!(zind,iz)
#         Ev̂ = Ev̂ + Π_z[:,iz] .* ( row_kron( Φ_z[ zind, :] , Φ_p̃ ) * cv)
#     end
#
#     val_adjust = profit(p̃, z_vals, w) + β * Ev̂
#
#     V[:,1] = val_adjust[ind_p_x_z[:,2]]
# end
#
# function value_nadjust!(V, grid_nodes, w, v̂)
#
#     @getPar __pars
#
#     #== Find the continuation expected value ==#
#     Ev̂ = v2Ev(v̂, Π_z)
#
#     p̃m1     = grid_nodes[:,1]
#     z_nodes = grid_nodes[:,2]
#
#     V[:,2] = profit(p̃m1, z_nodes, w) + β * Ev̂
#
#     return Void
# end
#
# function profit(p̃, z_nodes, w)
#     @getPar __pars
#
#     return ( exp(p̃) - w./ z_nodes ) .* ( exp( -ϵ * p̃)  )
# end
#
# function v2Ev(v,Π_z)
#
#     n_z = size(Π_z,1)
#     n_p = div(length(v),n_z )
#
#     return vec( reshape(v, n_p, n_z) * Π_z' )
#
# end
#
# function v2Ev!(Ev, v, Π_z)
#
#     n_z = size(Π_z,1)
#     n_p = div(length(v),n_z )
#
#     copy!(Ev, vec(reshape(v, n_p, n_z) * Π_z') )
#
# end
# # function profit_deriv(p̃, z_vals,w, Y, fcoll)
# #
# #     return ( exp(p̃) - w./ z_vals ) .* ( exp(p̃) .^( -ϵ ) ) * Y
# # end
# # --------------------------------------------------------------------------------------
# # % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# # --------------------------------------------------------------------------------------
