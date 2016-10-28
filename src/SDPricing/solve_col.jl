
# **************************************************************************************
#   REVIEW of CompEcon
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
# #       Φ = BasisStructure(basis), but doing it this way gives me the direct
# #       representation so I get Φ_y without repeating any calculations
# Φ_tensor = BasisStructure(basis,Tensor())
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
# Φ_direct = BasisStructure(basis,Direct(), S , 0)
#
# Φ_y  = Array( Φ_direct.vals[1] )
# Φ_y2 = Array( BasisStructure(y_basis, Expanded()).vals[1] )
# Φ = BasisStructure(basis,Expanded())
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
# Φ_new  = BasisStructure(basis, Expanded(), Snew, 0).vals[1]
# [Φ_new*coeff funeval(coeff, basis, Snew)]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function solve_firm_policy(w::Float64,Y::Float64,fp::FirmProblem)

    Φ_fac = fp.Φ_fac

    V = zeros(, 2)
    while
        bellman_rhs!(V, coeff, w, Y, fp)
    end

    coeff[:,1] = Φ_fac\V[:,1]
    coeff[:,2] = Φ_fac\V[:,2]
    coeff[:,3] = Φ_fac\V[:,3]
end

"""
Computes the RHS of our system of equations for a given guess of collocation
coefficients.

## Arguments

##### Basis Info
   - `basis::CompEcon.Basis`    : basis structure
   - `S`                        : collocation nodes (states + forecast consumption)
   - `coeff`                    : coefficient on basis function

##### Model
   - `model`        : Instance of Economy

## Returns
Return depend on the inputs... The basic function returns RHS evaluated at the collocation nodes only...

"""
function bellman_rhs!(V::Array{Float64,2}, coeff::Array{Float64,2}, w::Float64, Y::Float64, fp::FirmProblem)

    ca = coeff[:,1]
    cn = coeff[:,2]
    cv = coeff[:,3]

    #== Basis matrices ==#
    Φ       = fp.mbasis.Φ
    Φ_z,_   = fp.mbasis.Φ_tensor.vals

    z_vals   = fp.mbasis.z_nodes
    n_z     = length(z_vals)
    grid_nodes = fp.mbasis.grid_nodes

    #== Indices ==#
    ind_z_x_z = fp.mbasis.ind_z_x_z
    ind_z_x_p̃ = fp.mbasis.ind_z_x_p̃

    #==  ==#
    p̃_basis = fp.mbasis.p̃_basis

    #== FIND optimal p̃ in case of adjustment ==#
    f!(p̃, resid) = foc_price_adjust(resid, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
    g!(p̃, J)     = soc_price_adjust(resid, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
    @time res = nlsolve(f!, g!, ones(n_z), autodiff = true)
    p̃star = res.zero

    #== Value function ==#
    value_adjust!(V, p̃, z_vals, w, Y, cv, Π_z, p̃_basis, Φ_z, ind_z_x_z, ind_z_x_p̃)
    value_nadjust!(V, grid_nodes, w, Y, cv, Π_z, Φ)
end

function foc_price_adjust{T<:Real}(resid::Vector{T}, p̃::Vector{T}, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)

    n_z = length(z_vals)
    # .....................................................................................

    Φ_p̃_deriv = BasisStructure( p̃_basis, Direct(), exp(p̃), 1).vals[1]
    Φ = row_kron( Φ_p̃_deriv[ ind_z_x_z[:,2], :], Φ_z[ ind_z_x_z[:,1], :] )
    ∂v̂ = Φ * cv

    E∂v̂ = row_kron( eye(n_z) , fp.Π_z ) * ∂v̂

    for iz =1:n_z
        resid[iz] = exp( p̃[iz]) * Y - fp.ϵ * ( exp(p̃[iz]) - w/z_vals[iz] * Y +
        fp.β * exp( (fp.ϵ + 1.0) * p̃[iz]) * E∂v̂[iz]
    end

    return Void
end

function soc_price_adjust{T<:Real}(J::Matrix{T}, p̃::Vector{T}, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)

    n_z = length(z_vals)
    # .....................................................................................

    Φ_p̃_deriv  = BasisStructure( p̃_basis, Direct(), exp(p̃), 1).vals[1]
    Φ_p̃_deriv2 = BasisStructure( p̃_basis, Direct(), exp(p̃), 2).vals[1]
    Φ = row_kron( Φ_p̃_deriv[ ind_z_x_z[:,2], :], Φ_z[ ind_z_x_z[:,1], :] )
    ∂v̂ = Φ * cv

    Φ = row_kron( Φ_p̃_deriv2[ ind_z_x_z[:,2], :], Φ_z[ ind_z_x_z[:,1], :] )
    ∂v̂² = Φ * cv

    E∂v̂  = row_kron( eye(n_z) , fp.Π_z ) * ∂v̂
    E∂v̂² = row_kron( eye(n_z) , fp.Π_z ) * ∂v̂²

    for iz =1:n_z
        J[iz,iz] = exp( p̃[iz]) * Y - fp.ϵ * exp(p̃[iz]) * Y +
         fp.β * (fp.ϵ + 1.0) * exp( (fp.ϵ + 1.0) * p̃[iz]) * E∂v̂[iz] +
         fp.β * exp( (fp.ϵ + 2.0) * p̃[iz]) * E∂v̂²[iz]
    end

    return Void
end


function foc_price_adjust(p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)

    n_z = length(z_vals)
    # .....................................................................................

    Φ_p̃_deriv = BasisStructure( p̃_basis, Direct(), exp(p̃), 1).vals[1]
    Φ = row_kron( Φ_p̃_deriv[ ind_z_x_z[:,2], :], Φ_z[ ind_z_x_z[:,1], :] )
    ∂v̂ = Φ * cv

    E∂v̂ = row_kron( eye(n_z) , fp.Π_z ) * ∂v̂

    return exp(-fp.ϵ * p̃) * Y - fp.ϵ * ( exp(p̃) - w/z_vals ) * exp( -(fp.ϵ + 1.0) * p̃) * Y + fp.β * E∂v̂
end

function value_adjust!(V::Array{Float64,2}, p̃::Vector{Float64}, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z, ind_z_x_p̃)

    ##  NOTE:  Φ_z dimensions are nNodes × nAshock ##
    Φ_p̃ = BasisStructure( p̃_basis, Direct(), p̃, 0).vals[1]
    Φ   = row_kron( Φ_p̃[ ind_z_x_z[:,2], :] , Φ_z[ ind_z_x_z[:,1], :] )

    v̂ = Φ * cv
    Ev̂ = row_kron( eye(n_z), Π_z ) * v̂

    val_zdjust = profit(p̃, z_vals, w, Y) + fp.β * Ev̂

    V[:,1] = val_zdjust[ind_z_x_p̃[:,1]]
end

function value_nadjust!(V, grid_nodes, w, Y, cv, Π_z, Φ)

    n_z = size(Π_z,1)
    n_p̃ = div(size(grid_nodes,1), n_z)

    v̂ = Φ * cv
    Ev̂ = kron( eye(n_p̃) , Π_z ) * v̂

    p̃m1    = grid_nodes[:,1]
    z_vals  = grid_nodes[:,2]

    val_nadjust = profit(p̃m1, z_vals, w, Y) + fp.β * Ev̂

    V[:,2] = val_nadjust
end

function profit(z_vals, p̃, w, Y, fp)

    return ( p̃ - w./ z_vals ) .* ( p̃ .^( -fp.ϵ ) ) * Y
end
# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------
