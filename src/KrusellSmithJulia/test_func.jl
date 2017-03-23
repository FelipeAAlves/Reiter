# **************************************************************************************
#   TESTING policies
#
# ======================================================================================
#
# using NLsolve, Roots
#
# #== Initialize Objects ==#
# Cons         = Reiter.ConsumerProblem();
# ss_histogram = Reiter.StstHistogram(Cons);
#
# #== GUESS prices ==#
# Rss      = 1 + Reiter.netintr(1.01*Cons.KrepSS, 0);
# wagess   = Reiter.wagefunc(1.01*Cons.KrepSS, 0);
#
# ##  TEST eulerres fnc  ##
# out1 = zeros(120*2)
# out2 = zeros(120*2)
# @time Reiter.eulerres1!(out1, Cons.Θinit , Cons.Θinit, Rss, Rss, wagess, wagess, Cons)
# @time Reiter.eulerres2!(out2, Cons.Θinit , Cons.Θinit, Rss, Rss, wagess, wagess, Cons)
#
# @code_warntype Reiter.eulerres1!(out1, Cons.Θinit , Cons.Θinit, Rss, Rss, wagess, wagess, Cons)
# @code_warntype Reiter.eulerres2!(out2, Cons.Θinit , Cons.Θinit, Rss, Rss, wagess, wagess, Cons)
#
# ##  TEST policy function  ##
# f!(Θ, fvec) = Reiter.eulerres1!(fvec, Θ , Θ, Rss, Rss, wagess, wagess, Cons)
# g!(Θ, fvec) = Reiter.eulerres2!(fvec, Θ , Θ, Rss, Rss, wagess, wagess, Cons)
#
# #== Compare performance ==#
# # @time nlsolve(f!, Cons.Θinit)
# @time nlsolve(f!, Cons.Θinit, autodiff = true)
# @time nlsolve(g!, Cons.Θinit, autodiff = true)
# #
# res1 = nlsolve(f!, Cons.Θinit, autodiff = true)
# res2 = nlsolve(g!, Cons.Θinit, autodiff = true)
# # copy!(Cons.Θinit, res2.zero);

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# **************************************************************************************
#   TEST Derivatives
#
# ======================================================================================

# using ForwardDiff
#
# ## TEST 01 : test the jacobian on the INTEPOLATION ##
# out = zeros(120*2);
# f!(y, Θ) = Reiter.test_jaco_interp!(y, Θ, Rss, wagess, Cons)
# J = ForwardDiff.jacobian( f!, out, Cons.Θinit, Chunk{10}())
#
# @time ForwardDiff.jacobian( f!, out, Cons.Θinit, Chunk{2}() );
# @time ForwardDiff.jacobian( f!, out, Cons.Θinit, Chunk{5}() );
# @time ForwardDiff.jacobian( f!, out, Cons.Θinit, Chunk{10}());
#
# ##  TEST 02 : Do DUAL by hand  ##
# Θ     = Dual{2,Float64}[Dual(Cons.Θinit[i],i==1,i==2 ) for i=1:240]
# Θnext = Dual{2,Float64}[Dual(Cons.Θinit[i],0,0) for i=1:240]
#
# num, cons = Reiter.eulerres_test(Θ, Θnext, Rss, Rss, wagess, wagess, Cons)
# # CHECK the value with what comes out of jacobian function
# 1-num[2]/cons[2]
#
# f!(y, Θall) = Reiter.eulerres2!(y, Θall[1:240], Θall[241:end], Rss, Rss, wagess, wagess, Cons)
# @time J = ForwardDiff.jacobian(f!, out, repmat(Cons.Θinit,2) );
# # .....................................................................................
#
# ## TEST 03 : test the jacobian on Euler Residuals ##
# f!(y, Θall) = Reiter.eulerres2!(y, Θall[1:240], Θall[241:end], Rss, Rss, wagess, wagess, Cons)
#
# #== No pre-allocation ==#
# @time J = ForwardDiff.jacobian(f!, out, repmat(Cons.Θinit,2) );
#
# #== Pre-allocated ==#
# Jout = ForwardDiff.JacobianResult(out);
# @time ForwardDiff.jacobian!(Jout, f!, Cons.Θinit, Chunk{10}());
# @time ForwardDiff.jacobian!(Jout, f!, Cons.Θinit, Chunk{5}());
# @time ForwardDiff.jacobian!(Jout, f!, Cons.Θinit, Chunk{2}());
#
# Jout = similar(J)
# f(Θ) = eulerres(Cons.Θinit, Θ, Rss, Rss, wagess, wagess, Cons)
# @time J = ForwardDiff.jacobian!(Jout, f, Cons.Θinit);
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# **************************************************************************************
#   TEST Histogram SteadyState
# ======================================================================================
using Roots

#== Initialize Objects ==#
Cons         = Reiter.ConsumerProblem();
ss_histogram = Reiter.StstHistogram(Cons);

f_hist(k) = Reiter.stst_histogram_resid(Cons, k)

#== ProfileView ==#
# f_hist(1.01*Cons.KrepSS)
# Profile.clear()
# @profile f_hist(1.01*Cons.KrepSS)
# using ProfileView
# ProfileView.view()

Kss = fzero(f_hist, 1.01*Cons.KrepSS, 1.05*Cons.KrepSS)

#== UPDATE ss_histogram ==#
Reiter.stst_histogram_resid(Cons, Kss, ss_histogram)

# **************************************************************************************
#   TEST Density SteadyState
# ======================================================================================
# ss_density = Reiter.StstDensity(Cons);
#
# #== Compute moments from histogram  ==#
# mMoments = zeros(Cons.nMoments, Cons.nEps);
# for ieps = 1:Cons.nEps
#     for iMom = 1:Cons.nMoments
#         dens = ss_histogram.mHistogram[:,ieps]/sum( ss_histogram.mHistogram[:,ieps] )
#         mMoments[iMom, ieps] = dot( (Cons.asset_grid_fine - (iMom>1) * mMoments[1, ieps] ).^iMom, dens )
#     end
# end
#
# f_dens(k) = Reiter.stst_density_resid(Cons, k, mMoments)
# Kss_density = fzero(f_dens, 0.995*ss_histogram.capital, 1.005*ss_histogram.capital)
#
# #== UPDATE ss_density ==#
# Reiter.stst_density_resid(Cons, Kss_density, mMoments, ss_density)

# ******************************************************************
#   FIGURE : Histogram × Density
# ==================================================================
# using Plots, LaTeXStrings
# pyplot() # plotlyjs()
#
# asset_grid_fine = Cons.asset_grid_fine;
# asset_quad_points = Cons.asset_quad_points;
# Δa = asset_grid_fine[2]-asset_grid_fine[1];
# z_dist = Cons.z_dist
#
# #== CHECK density ==#
# Reiter.density_int(ss_density.mDensityCoeff[2:end,1], ss_density.mMoments[:,1], Cons, 0, ss_density.mDensityCoeff[1,1])
# Reiter.density_int(ss_density.mDensityCoeff[2:end,2], ss_density.mMoments[:,2], Cons, 0, ss_density.mDensityCoeff[1,2])
#
# ieps = 1
#     p1 = bar(asset_grid_fine, ss_histogram.mHistogram[:,ieps]/(Δa*z_dist[ieps]), alpha=0.6, label=LaTeXString("\$ h_1(a_i)\$"),  legendfont=font(12));
#     dens = ss_density.mDensityCoeff[1,ieps] * Reiter.density(ss_density.mDensityCoeff[2:end,ieps], ss_density.mMoments[:,ieps], asset_grid_fine, Cons);
#     plot!(asset_grid_fine, dens, color=:red, linewidth=2, label=LaTeXString("\$ g_{1} (a) \$"),  legendfont=font(12));
# ieps = 2
#     p2 = bar(asset_grid_fine, ss_histogram.mHistogram[:,ieps]/(Δa*z_dist[ieps]), alpha=0.6,label=LaTeXString("\$ h_2(a_i, 1_j)\$"),  legendfont=font(12));
#     dens = ss_density.mDensityCoeff[1,ieps] * Reiter.density(ss_density.mDensityCoeff[2:end,ieps], ss_density.mMoments[:,ieps], asset_grid_fine, Cons);
#     plot!(asset_grid_fine, dens, color=:red, linewidth=2, label=LaTeXString("\$g_{2}(a)\$"), legendfont=font(12));
#
# titles = ["Unemp" "Employed"]
# p_hist_x_dens = plot(p1, p2, layout=2, title=titles);

# **************************************************************************************
#   Linearization STEP
# ======================================================================================
println("STEP [2]: Linearizing model equations")


#=======================================================#
###          CONSTRUCT iZvar for use in the          ###
#=======================================================#

iZvar = Dict{Symbol, UnitRange{Int64}}()

nExog = 1
nx,ny = length(ss_histogram.xstst), length(ss_histogram.ystst)

#== CREATE Y ==#
yxstst = [ss_histogram.ystst; ss_histogram.xstst];

Zss    = [ yxstst; yxstst; zeros(nExog) ];

#== iZvar ==#
iZvar[:y′]    = 1:ny
iZvar[:x′]    = ny+1:ny+nx

iZvar[:y] = (nx+ny) + (iZvar[:y′])
iZvar[:x] = (nx+ny) + (iZvar[:x′])
nn = 2*(nx+ny)
iZvar[:eps] = nn+1:nn+nExog


#== Evaluate Stst ==#
# using Gallium
# Gallium.breakpoint(Reiter.equil_histogram)

# equil_res = Reiter.equil_histogram(Yss, iYvar, ss_histogram, Cons)
equil_res = Reiter.equil_histogram(Zss, iZvar, ss_histogram, Cons)
maximum( abs( equil_res) ) > 1e-6  && throw(error("Wrong Stst"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
using ForwardDiff
using Gensys
using QuantEcon

# f(Y) = Reiter.equil_histogram(Y, iYvar, ss_histogram, Cons)
f(Z) = Reiter.equil_histogram(Z, iZvar, ss_histogram, Cons)

∇f = zeros( length(equil_res), length(Zss) );
@time ForwardDiff.jacobian!(∇f, f, Zss);

∇f_y′     =   ∇f[:,  iZvar[:y′] ] ;
∇f_x′     =   ∇f[:,  iZvar[:x′] ] ;
∇f_y      =   ∇f[ :, iZvar[:y] ]  ;
∇f_x      =   ∇f[ :, iZvar[:x] ]  ;
Ψ         = - ∇f[ :, iZvar[:eps] ];
η_shocks  =  Ψ[ ny+1:end, : ] # equations refering to transition of states

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# using Gallium
# Gallium.breakpoint(Reiter.klein)

#============================#
##           KLEIN          ##
#============================#
# gx, hx = Reiter.klein(∇f_x′, ∇f_y′, ∇f_x, ∇f_y, 1.0+1e-10)

#== Adjust the exogenous shock ==#
# ∇f_x[:,ss_histogram.ixvar[:aggr_exo]] = ∇f_x[:,ss_histogram.ixvar[:aggr_exo]] + Reiter.ρz*∇f_x′[:,ss_histogram.ixvar[:aggr_exo]]
# ∇f_x′[:,ss_histogram.ixvar[:aggr_exo]] = 0.0

_nVAR = nx + ny;
    Γ0 =  [ ∇f_x′ ∇f_y′ ];
    Γ1 = -[ ∇f_x ∇f_y ];
    c  = zeros(_nVAR);


# writedlm("gamma0.txt", Γ0)
# Γ0_ = readdlm("gamma0.txt")

# using Gallium
# Gallium.breakpoint(gensysdt)
# Basic Commands:
#
# n steps to the next line
# s steps into the next call
# finish runs to the end of the function
# bt shows a simple backtrace
# `stuff runs stuff in the current frame's context
# fr v will show all variables in the current frame
# f n where n is an integer, will go to the n-th frame.

println("STEP [3]: LINEARIZING the model")
@time Phi1, cons, B1, _, _, _, _, eu1, _ = gensysdt(Γ0, Γ1, c , Ψ , Π1, 1.0 + 1e-10);
@time Phi2, cons, B2, _, _, _, _, eu2, _ = gensysdt(Γ0, Γ1, c , Ψ , Π2, 1.0 + 1e-10);
eu1
eu2

nx,ny = length(ss_histogram.xstst), length(ss_histogram.ystst)
ϵ_sim = randn(1,10)

x_sim1 = zeros(_nVAR,11)
x_sim2 = zeros(_nVAR,11)
for it =1:10
    x_sim1[1:nx,it+1] = hx * x_sim1[1:nx,it] + η_shocks*ϵ_sim[:,it];
    x_sim1[nx+1:end,it+1] = gx * x_sim1[1:nx,it+1];
    x_sim2[:,it+1] = Phi2*x_sim2[:,it] + B2 * ϵ_sim[:,it];
end
maxabs(x_sim1[:] - x_sim2[:])


# **************************************************************************************
#   FIGURES
#
# ======================================================================================
# """
# Linear State Space Models
#     x_{t+1} = Phi x_t + B w_{t+1}
#         y_t = H x_t
#
# - `Phi::Matrix` Part of the state transition equation.  It should be `n x n`
# - `B::Matrix` Part of the state transition equation.  It should be `n x m`
# - `H::Matrix` Part of the observation equation.  It should be `k x n`
# - `k::Int` Dimension
# - `n::Int` Dimension
# - `m::Int` Dimension
# - `mu_0::Vector` This is the mean of initial draw and is of length `n`
# - `Sigma_0::Matrix` This is the variance of the initial draw and is `n x n` and
#                     also should be positive definite and symmetric
#
# """

# H   = eye(size(Phi,1));
# lss = LSS(Phi, B, H)  ;# mu_0 = zeros(size(G, 2), Sigma_0 = zeros(size(H, 2), size(H, 2)) )
#
# simul_length = 200
# X_simul, _ = simulate(lss, simul_length);
# X_simul = X_simul + repmat(Xss,1,simul_length);
#
# nEps        = Cons.nEps
# nSavingsPar = Cons.nSavingsPar
# nAssetsFine = Cons.nAssetsFine
# asset_grid_fine = Cons.asset_grid_fine
# nHistogram  = nAssetsFine * Cons.nEps
#
# # WARN : Change
# iXvar              = ss_histogram.iXvar
#
# K_simul            = X_simul[ iXvar[:aggr_end],: ]';
# z_simul            = X_simul[ iXvar[:aggr_exo],: ]';
#
# Θ_simul            = X_simul[ iXvar[:household],: ];
# Cpol_simul         = zeros(nAssetsFine*nEps, simul_length);
# HistogramDev_simul = X_simul[ iXvar[:histogram],: ];
# Histogram_simul    = zeros(nAssetsFine,simul_length);
#
# vHistogram  = zeros(nHistogram) ;
# vHistogramK = zeros(nAssetsFine);
# Cpol        = zeros(nAssetsFine*nEps);
#
# for t = 1:simul_length
#
#     #== Histogram ==#
#     Reiter.x2distr!(vHistogram, HistogramDev_simul[:,t], ss_histogram.vHistogram)
#     copy!(vHistogramK, sum( reshape(vHistogram, nAssetsFine, 2) , 2) )
#     Histogram_simul[:,t] = vHistogramK
#
#     #== Consumption ==#
#     R      = 1 + Reiter.netintr(K_simul[t], z_simul[t]);
#     wage   = Reiter.wagefunc(K_simul[t], z_simul[t]);
#     for ieps = 1:nEps
#         ind = (ieps-1)*nSavingsPar+1 : ieps*nSavingsPar
#         Sthis= Reiter.interpolate(Cons, Θ_simul[ind,t])
#         sthis = Sthis[asset_grid_fine]
#
#         ind = (ieps-1)*nAssetsFine+1 : ieps*nAssetsFine
#         Cpol[ind] = Reiter.get_consumption(asset_grid_fine, sthis, R, wage, ieps, Cons)
#     end
#     Cpol_simul[:,t] = Cpol
# end
#
# #== HIGH and Low period ==#
# ind = findmax(z_simul)[2], findmin(z_simul)[2];
#
# Cpol_SS = zeros(nAssetsFine, nEps)
# for ieps =1:nEps
#     Sthis = Reiter.interpolate(Cons, ss_histogram.mΘ[:,ieps])
#     sthis = Sthis[asset_grid_fine]
#     Cpol_SS[:,ieps] = Reiter.get_consumption(asset_grid_fine, sthis, ss_histogram.R, ss_histogram.wage, ieps, Cons)
# end
#
# # .....................................................................................
# plotlyjs()
# y_vals = []
# labels = []
#
# push!(y_vals, Cpol_SS[:,1])
# push!(y_vals, Cpol_simul[1:nAssetsFine, ind[1] ])
# push!(y_vals, Cpol_simul[1:nAssetsFine, ind[2] ])
#
# push!(y_vals, Cpol_SS[:,2])
# push!(y_vals, Cpol_simul[(1:nAssetsFine)+nAssetsFine, ind[1] ])
# push!(y_vals, Cpol_simul[(1:nAssetsFine)+nAssetsFine, ind[2] ])
#
# #== Consumption policies ==#
# psim1 = plot(asset_grid_fine, y_vals, linewidth=2, alpha=0.6, xlim=[0.0, 0.5*asset_grid_fine[end] ]);
#
# #== Histogram ==#
# psim2 = bar(Cons.asset_grid_fine, [Histogram_simul[:,1] Histogram_simul[:,ind[1]] Histogram_simul[:,ind[2]]], alpha=0.6);
#
# # commands to plot
# # plot(psim1)
# # plot(psim2)
# # plot(p_hist_x_dens)
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
