
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# **************************************************************************************
#   FIND optimal p̃ in case of adjustment
#
# ======================================================================================
using BasisMatrices
using NLsolve

# fcoll = Reiter.FirmColloc(Spline,Spline)
# fpol = Reiter.FirmPolicies(fcoll)
# fcoll = Reiter.FirmColloc(Spline,Cubic)

coeff = fcoll.coeff;
ca = coeff[:,1];
cn = coeff[:,2];
cv = coeff[:,3];

#== Basis matrices ==#
Φ       = fcoll.Φ;
Φ_z,_   = fcoll.Φ_tensor.vals;
grid_nodes = fcoll.grid_nodes;

#== Indices ==#
ind_z_x_z = fcoll.ind_z_x_z;
ind_p_x_z = fcoll.ind_p_x_z;

#==  ==#
p_basis = fcoll.p_basis;

#== Guesses ==#
w = 1.0

Reiter.@getPar(Reiter.__pars)

pstar = log( ϵ/(ϵ-1.0) * w./z_vals);
resid = zeros(10);
Reiter.foc_price_adjust!(resid, zeros(10), w, cv, p_basis, Φ_z, ind_z_x_z)
resid
Reiter.foc_price_adjust!(resid, pstar, w, cv, p_basis, Φ_z, ind_z_x_z)
resid
Reiter.foc_price_adjust!(resid, fpol.pstar, w, cv, p_basis, Φ_z, ind_z_x_z)
resid

f!(p̃, fvec) = Reiter.foc_price_adjust!(fvec, p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
g!(p̃, J)    = Reiter.soc_price_adjust!(J, p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
f_val(p̃) = Reiter.foc_price_adjust(p̃, w, cv, p_basis, Φ_z, ind_z_x_z)
f_jac(p̃) = Reiter.soc_price_adjust(p̃, w, cv, p_basis, Φ_z, ind_z_x_z)

@time (pstar3, x) = Reiter.broydn(f_val, zeros(10),[1e-6,1,2],f_jac)
@time res1 = nlsolve(f!,  zeros(10))
@time res2 = nlsolve(f!, g!,  zeros(10))
[exp(res1.zero) exp(res2.zero) exp(pstar3) exp(fpol.pstar)]

# #== FIND optimal p̃ in case of adjustment ==#
#
# # Solve for optimal price
# pstar = log( fcoll.ϵ/(fcoll.ϵ-1.0) * w./z_vals);
#
# resid = zeros(10);
# Reiter.foc_price_adjust!(resid, zeros(10), z_vals, w, Y, cv, fcoll, p_basis, Φ_z, ind_z_x_z)
# resid
#
# f!(p̃, fvec) = Reiter.foc_price_adjust!(fvec, p̃, z_vals, w, Y, cv, fcoll, p_basis, Φ_z, ind_z_x_z)
# @time res = nlsolve(f!, zeros(10))
# [exp(res.zero) exp(pstar)]
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# # test function
# resid = zeros(10);
# Reiter.foc_price_adjust!(resid, zeros(10), z_vals, w, Y, cv, fcoll, p_basis, Φ_z, ind_z_x_z)
# resid
#
# resid[:] = 0.0;
# Reiter.foc_price_adjust!(resid, pstar, z_vals, w, Y, cv, fcoll, p_basis, Φ_z, ind_z_x_z)
# resid
#
# Jout = zeros(10,10);
# Reiter.soc_price_adjust!(Jout, pstar, z_vals, w, Y, cv, fcoll, p_basis, Φ_z, ind_z_x_z)
# Jout
#
# f!(p̃, fvec) = Reiter.foc_price_adjust!(fvec, p̃, z_vals, w, Y, cv, fcoll, p_basis, Φ_z, ind_z_x_z)
# g!(p̃, J)    = Reiter.soc_price_adjust!(J, p̃, z_vals, w, Y, cv, fcoll, p_basis, Φ_z, ind_z_x_z)
#
# @time res = nlsolve(f!, zeros(10), autodiff = true)
# @time res1 = nlsolve(f!, zeros(10))
# @time res2 = nlsolve(f!, g!, zeros(10))
# [exp(res1.zero) exp(res2.zero) exp(pstar)]


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------
using BasisMatrices

fcoll_1 = Reiter.FirmColloc(Spline,Spline)
sol_1   = Reiter.FirmSolution(fcoll_1)

fcoll_2 = Reiter.FirmColloc(Spline,Spline)
sol_2   = Reiter.FirmSolution(fcoll_2)
w = 0.5
Y = 1.0

Φ = fcoll.Φ
V    = Φ * sol.coeff
V1 = similar(V)
V2 = similar(V)

@time Reiter.bellman_rhs!(V1, fcoll, sol, w, Y)
@time Reiter.bellman_rhs_eval_v!(V2, fcoll, sol, w, Y)


@time Reiter.solve_firm_policy!(fcoll_1, sol_1, w, Y, 1)
@time Reiter.solve_firm_policy!(fcoll_2, sol_2, w, Y, 2)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Reiter.@getPar(Reiter.__pars)

#== Basis matrices ==#
Φ       = fcoll.Φ
Φ_z     = fcoll.Φ_tensor.vals[2]
p_basis = fcoll.p_basis

#== Indices ==#
grid_nodes  = fcoll.grid_nodes
ind_z_x_z   = fcoll.ind_z_x_z
ind_p_x_z   = fcoll.ind_p_x_z

#== Coefficients for before shock value function ==#
cv = sol.coeff[:,3]
v̂  = Φ * cv            # continuation value

V    = Φ * sol.coeff
V1    = similar(V)
V2    = similar(V)
pstar = log( ϵ/(ϵ-1) * w ./ z_vals )

function eval_v̂(p̃, ind_a)

    @assert length(p̃)==length(ind_a)
    Φ_p̃ = BasisMatrix( p_basis, Direct(), p̃, 0).vals[1]
    Φ   = row_kron( Φ_z[ ind_a, :],  Φ_p̃ )
    return Φ * cv
end
@time Reiter.value_adjust!(V1, pstar, w, Y, cv, p_basis, Φ_z, ind_z_x_z, ind_p_x_z);
@time Reiter.value_adjust_test!(V2, pstar, eval_v̂,  w, Y, ind_z_x_z, ind_p_x_z);


# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

# **************************************************************************************
#   USING code
#
# ======================================================================================
using BasisMatrices
using NLsolve
fcoll = Reiter.FirmColloc(Spline)
sol   = Reiter.FirmSolution(fcoll)

# fill!(sol.coeff,0.0);
# fill!(sol.pstar,0.0);
# Reiter.solve_firm_policy!(1.0, fcoll, sol)
#
# fill!(sol.coeff,0.0);
# fill!(sol.pstar,0.0);
# Reiter.solve_firm_policy!(0.5, fcoll, sol)
#
# fill!(sol.coeff,0.0);
# fill!(sol.pstar,0.0);
# Reiter.solve_firm_policy!(0.30, fcoll,sol)


#==========================#
##  Finding the SS        ##
#==========================#
# using Roots

ss_histogram = Reiter.StstHistogram(fcoll)

##  WARN : check why result is not changing  ##
# Reiter.stst_histogram_resid(0.75, 1.00, ss_histogram, fcoll, sol)
# Reiter.stst_histogram_resid(1.25, 1.00,ss_histogram, fcoll, sol)


function f_hist!(x,fvec)
    copy!( fvec, Reiter.stst_histogram_resid(x[1], x[2], ss_histogram, fcoll, sol) )
end

res = nlsolve(f_hist!, [0.75,1.0/3], iterations = 50)
Reiter.stst_histogram_resid(res.zero[1], res.zero[2], ss_histogram, fcoll, sol, true)

# reshape(ss_histogram.Xss[ss_histogram.iXvar[:value_fnc]],32,10)

#=============================================#
###   CONSTRUCT iZvar for use in res_equi   ###
#=============================================#
iZvar = Dict{Symbol, UnitRange{Int64}}()

nExog = 2
# nEta  = length(fcoll.grid_nodes[:,1]) + length(ss_histogram.vHistogram) + n_z
nx,ny = length(ss_histogram.xstst), length(ss_histogram.ystst)

#== CREATE Y ==#
yxstst = [ss_histogram.ystst; ss_histogram.xstst];
Zss = [ yxstst; yxstst; zeros(nExog) ];

#== iZvar ==#
iZvar[:y′]    = 1:ny
iZvar[:x′]    = ny+1:ny+nx

iZvar[:y] = (nx+ny) + (iZvar[:y′])
iZvar[:x] = (nx+ny) + (iZvar[:x′])
nn = 2*(nx+ny)

iZvar[:eps] = nn+1:nn+nExog
# nn +=nExog
# iZvar[:eta] = nn+1:nn+nExog

# .....................................................................................

# using Gallium
# Gallium.breakpoint(Reiter.equil_histogram)

res_equil = Reiter.equil_histogram(Zss, iZvar, ss_histogram, fcoll, sol);
maximum( abs( res_equil) ) > 1e-6  && throw(error("Wrong Stst"))
_nVAR = length(res_equil)

# @code_warntype Reiter.equil_histogram(Zss, iZvar, ss_histogram, fcoll, sol)

# using Gallium
# Gallium.breakpoint(Reiter.equil_histogram)

# Basic Commands:
#
# n steps to the next line
# s steps into the next call
# finish runs to the end of the function
# bt shows a simple backtrace
# `stuff runs stuff in the current frame's context
# fr v will show all variables in the current frame
# f n where n is an integer, will go to the n-th frame.



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
using ForwardDiff

f(Y) = Reiter.equil_histogram(Y, iZvar, ss_histogram, fcoll, sol)
∇f = zeros( length(res_equil), length(Zss) );
@time ForwardDiff.jacobian!(∇f, f, Zss, ForwardDiff.JacobianConfig{5}(Zss));
# @time ∇f = ForwardDiff.jacobian(f, Zss, ForwardDiff.JacobianConfig{5}(Zss))


#============================#
###   Applying Klein algorithm   ###
#============================#

function klein(∇f_x′, ∇f_y′, ∇f_x, ∇f_y, stake=1.0)

    Γ0 = -[ ∇f_x′ ∇f_y′ ];
    Γ1 =  [ ∇f_x ∇f_y ];
    n_k = size(∇f_x,2);
    F = schurfact(complex(Γ0), complex(Γ1));
    #  Γ0=F[:Q]*F[:S]*F[:Z]' and Γ_1=F[:Q]*F[:T]*F[:Z]'
    klein(F, n_k, stake)
end

function klein(F::Base.LinAlg.GeneralizedSchur, n_k::Int64, stake::Float64)

    S, T = F[:S], F[:T]
    movelast = Bool[abs(T[i, i]) > stake * abs(S[i, i]) for i in 1:size(S,1)]
    n_stable = sum(!movelast)

    @printf("   Number of state variables    : %d \n", n_k)
    @printf("   Number of stable eigenvalues : %d \n", n_stable)
    if n_stable>n_k
        error("The Equilibrium is Locally Indeterminate")
    elseif n_stable<n_k
        error("No Local Equilibrium Exists")
    else
        FS = ordschur!(F, !movelast)
        S, T, Qt, Z = FS[:S], FS[:T], FS[:Q], FS[:Z]
        # [abs(T[i, i]/S[i, i])  for i in 1:size(S,1)]
        Z_11 = Z[1:n_stable,1:n_stable];
        Z_21 = Z[n_stable+1:end,1:n_stable];

        S_11 = S[1:n_stable,1:n_stable];
        T_11 = T[1:n_stable,1:n_stable];

        if rank(Z_11)<n_stable;
            error("Invertibility condition violated")
        else
            @printf("   Invertibility CHECKED\n")
        end

        Z_11i = inv(Z_11);
        gx = real(Z_21*Z_11i);
        hx = real(Z_11 * (S_11\T_11) * Z_11i );

        return gx, hx
    end
end

∇f_y′   = ∇f[:,  iZvar[:y′] ] ;
∇f_x′   = ∇f[:,  iZvar[:x′] ] ;
∇f_y    =  ∇f[ :, iZvar[:y] ]  ;
∇f_x    =  ∇f[ :, iZvar[:x] ]  ;
η_shocks    =  -∇f[ ny+1:end, iZvar[:eps] ]; # equations refering to transition of states

gx, hx = klein(∇f_x′, ∇f_y′, ∇f_x, ∇f_y)

# **************************************************************************************

Γ0 = -[ ∇f_x′ ∇f_y′ ];
Γ1 =  [ ∇f_x ∇f_y ];
c  = zeros(_nVAR);
Ψ  = ∇f[ :, iZvar[:eps] ];

_Pi_    = -∇f_y′; # E_t{ y_{t+1} } = y_t+1 - eta_t+1
ind_η  = Reiter.findnzcols(_Pi_, 1e-12) # Reiter.findnzcols(_Pi_);
_nPi = length(ind_η)
Π = _Pi_[:, ind_η ]
#== effective/contaminated foward looking equations ==#
# ind_Pi = Reiter.findnzrows(_Pi_)
# ind_Pi_effec = ind_Pi[1:_nPi]
# ind_Pi_extra = ind_Pi[_nPi+1:end]
#
# Λ = _Pi_[ind_Pi_effec,:]
# _Pi = spzeros(_nVAR,_nPi)
# _Pi[ind_Pi_effec,:] = speye(_nPi)
# _Pi[ind_Pi_extra,:] = _Pi_[ind_Pi_extra,ind_η]/Λ

using Gensys
println("STEP [3]: LINEARIZING the model")
@time Phi, cons, B, _, _, _, _, eu, _ = gensysdt(Γ0, Γ1, c , Ψ , Π, 1.0 + 1e-10);
eu


nx,ny = length(ss_histogram.xstst), length(ss_histogram.ystst)
ϵ_sim = randn(2,10)

x_sim1 = zeros(_nVAR,11)
x_sim2 = zeros(_nVAR,11)
for it =1:10
    x_sim1[1:nx,it+1] = hx * x_sim1[1:nx,it] + η_shocks*ϵ_sim[:,it];
    x_sim1[nx+1:end,it+1] = gx * x_sim1[1:nx,it+1];
    x_sim2[:,it+1] = Phi*x_sim2[:,it] + B * ϵ_sim[:,it];
end
maxabs(x_sim1[:] - x_sim2[:])
