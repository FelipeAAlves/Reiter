

# **************************************************************************************
#   TEST MAIN
#
#   - compare stst histogram × density
#   - simlulate the model for the steady state
#   - check the distribution moves around
# ======================================================================================
# import Reiter
using Reiter
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

#== Find Stst for histogram ==#
Kss = fzero(f_hist, 1.01*Cons.KrepSS, 1.05*Cons.KrepSS)
#== UPDATE ss_histogram ==#
Reiter.stst_histogram_resid(Cons, Kss, ss_histogram)

# **************************************************************************************
#   TEST Density SteadyState
# ======================================================================================
#== Initialize StstDensity ==#
ss_density = Reiter.StstDensity(Cons);

#== Compute moments from histogram  ==#
mMoments = zeros(Cons.nMoments, Cons.nEps);
for ieps = 1:Cons.nEps
    for iMom = 1:Cons.nMoments
        dens = ss_histogram.mHistogram[:,ieps]/sum( ss_histogram.mHistogram[:,ieps] )
        mMoments[iMom, ieps] = dot( (Cons.asset_grid_fine - (iMom>1) * mMoments[1, ieps] ).^iMom, dens )
    end
end

#== Find Stst for density ==#
f_dens(k) = Reiter.stst_density_resid(Cons, k, mMoments)
Kss_density = fzero(f_dens, 0.995*ss_histogram.capital, 1.005*ss_histogram.capital)

#== UPDATE ss_density ==#
Reiter.stst_density_resid(Cons, Kss_density, mMoments, ss_density)

# ******************************************************************
#   FIGURE : Histogram × Density
# ==================================================================
# Plot test
# using PyPlot
# x = linspace(0, 10, 200)
# y = sin(x)
# fig, ax = subplots()
# ax[:plot](x, y, "r-", linewidth=2, label=L"$y = \sin(x)$", alpha=0.6)
# ax[:legend](loc="upper center")

using Plots, LaTeXStrings
plotlyjs()

x = linspace(0, 10, 200)
y = sin(x)
Plots.plot(x, y,  label=string(L"$ h_1(a_i)$"), lw=2, alpha=0.6)
Plots.plot!(xlabel="current assets", ylabel="next period assets", grid=false)


asset_grid_fine = Cons.asset_grid_fine;
asset_quad_points = Cons.asset_quad_points;
Δa = asset_grid_fine[2]-asset_grid_fine[1];
z_dist = Cons.z_dist

#== CHECK density ==#
Reiter.density_int(ss_density.mDensityCoeff[2:end,1], ss_density.mMoments[:,1], Cons, 0, ss_density.mDensityCoeff[1,1])
Reiter.density_int(ss_density.mDensityCoeff[2:end,2], ss_density.mMoments[:,2], Cons, 0, ss_density.mDensityCoeff[1,2])

# #== with latex ==#
# using Plots, LaTeXStrings
# pyplot() #
# ieps = 1
#     p1 = bar(asset_grid_fine, ss_histogram.mHistogram[:,ieps]/(Δa*z_dist[ieps]), alpha=0.6, label=LaTeXString("\$ h_1(a_i)\$"),  legendfont=font(12));
#     dens = ss_density.mDensityCoeff[1,ieps] * Reiter.density(ss_density.mDensityCoeff[2:end,ieps], ss_density.mMoments[:,ieps], asset_grid_fine, Cons);
#     plot!(asset_grid_fine, dens, color=:red, linewidth=2, label=LaTeXString("\$ g_{1} (a) \$"),  legendfont=font(12));
# ieps = 2
#     p2 = bar(asset_grid_fine, ss_histogram.mHistogram[:,ieps]/(Δa*z_dist[ieps]), alpha=0.6,label=LaTeXString("\$ h_2(a_i, 1_j)\$"),  legendfont=font(12));
#     dens = ss_density.mDensityCoeff[1,ieps] * Reiter.density(ss_density.mDensityCoeff[2:end,ieps], ss_density.mMoments[:,ieps], asset_grid_fine, Cons);
#     plot!(asset_grid_fine, dens, color=:red, linewidth=2, label=LaTeXString("\$g_{2}(a)\$"), legendfont=font(12));
# titles = ["Unemp" "Employed"]
# p_hist_x_dens = plot(p1, p2, layout=2, title=titles);

#== using Plots ==#
###  NOTE:  use string to be able to add latex to it  ###
ieps = 1
    p1 = Plots.bar(asset_grid_fine, ss_histogram.mHistogram[:,ieps]/(Δa*z_dist[ieps]), alpha=0.6, label=string(L"$h_1(a_i)$"),  legendfont=font(12));
    dens = ss_density.mDensityCoeff[1,ieps] * Reiter.density(ss_density.mDensityCoeff[2:end,ieps], ss_density.mMoments[:,ieps], asset_grid_fine, Cons);
    plot!(asset_grid_fine, dens, color=:red, linewidth=2, label=string(L"$ g_{1} (a)$"),  legendfont=font(12));
ieps = 2
    p2 = Plots.bar(asset_grid_fine, ss_histogram.mHistogram[:,ieps]/(Δa*z_dist[ieps]), alpha=0.6,label=string(L"$ h_2(a_i) $"),  legendfont=font(12));
    dens = ss_density.mDensityCoeff[1,ieps] * Reiter.density(ss_density.mDensityCoeff[2:end,ieps], ss_density.mMoments[:,ieps], asset_grid_fine, Cons);
    plot!(asset_grid_fine, dens, color=:red, linewidth=2, label=string(L"$ g_{2}(a) $"), legendfont=font(12));

titles = ["Unemp" "Employed"]
p_hist_x_dens = Plots.plot(p1, p2, layout=2, title=titles);

# **************************************************************************************
#   LINEARIZATION STEP
#
#    - only for the histogram
# ======================================================================================

#== CHECK the steady-state ==#
equil_res = Reiter.equil_histogram(ss_histogram.Zstst, ss_histogram, Cons)
maximum( abs( equil_res) ) > 1e-6  && throw(error("Wrong Stst"))
_nVAR = length(equil_res);

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
using ForwardDiff

#== Compute the Jacobian ==#
f(Z) = Reiter.equil_histogram(Z, ss_histogram, Cons)
∇f = zeros( length(equil_res), length(ss_histogram.Zstst) );
@time ForwardDiff.jacobian!(∇f, f, ss_histogram.Zstst);

∇f_x′     =   ∇f[:,  ss_histogram.iZvar[:x′] ] ;
∇f_y′     =   ∇f[:,  ss_histogram.iZvar[:y′] ] ;
∇f_x      =   ∇f[ :, ss_histogram.iZvar[:x] ]  ;
∇f_y      =   ∇f[ :, ss_histogram.iZvar[:y] ]  ;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#--------------------------------#
##             KLEIN             ##
#--------------------------------#

    #== CHECK stable/unstable eigenvalues ==#
    Γ0 =  [  ∇f_x′ ∇f_y′];
    Γ1 = -[ ∇f_x ∇f_y ];
    writecsv("gamma0.csv", Γ0)
    writecsv("gamma1.csv", Γ1)


    n_k = size(∇f_x,2);
    F = schurfact(complex(Γ0), complex(Γ1));
    S, T = F[:S], F[:T];
    movelast = Bool[abs(T[i, i]) > 1.00 * abs(S[i, i]) for i in 1:size(S,1)]
    n_stable = sum(!movelast)

gx, hx = Reiter.klein(∇f_x′, ∇f_y′, ∇f_x, ∇f_y);

#----------------------------------#
##              SIMS              ##
#----------------------------------#

nx, ny = size(∇f_x,2),size(∇f_y,2)
#== Adjust the exogenous shock ==#
# ∇f_x[:,ss_histogram.ixvar[:aggr_exo]] = ∇f_x[:,ss_histogram.ixvar[:aggr_exo]] + Reiter.ρz*∇f_x′[:,ss_histogram.ixvar[:aggr_exo]]
# ∇f_x′[:,ss_histogram.ixvar[:aggr_exo]] = 0.0

    Γ0 =  [∇f_x′ ∇f_y′];
    Γ1 = -[∇f_x ∇f_y];
    Ψ = - ∇f[ :, ss_histogram.iZvar[:eps] ];
    c  = zeros(_nVAR);

    _Pi_  = ∇f_y′; # notation E_t{ y_{t+1} } = y_t+1 - eta_t+1
    ind_η = Reiter.findnzcols(_Pi_, 1e-12) # Reiter.findnzcols(_Pi_);
    _nPi  = length(ind_η)
    Π1    = _Pi_[:, ind_η ];
    Π2 = [eye(ny);zeros(nx,ny)];

@time Phi1, cons, B1, _, _, _, _, eu1, _ = Reiter.gensysdt(Γ0, Γ1, c , Ψ , Π1, 1.0 + 1e-10);
@time Phi2, cons, B2, _, _, _, _, eu2, _ = Reiter.gensysdt(Γ0, Γ1, c , Ψ , Π2, 1.0 + 1e-10);
eu1
eu2


# **************************************************************************************
#   SIMULATION
#
# ======================================================================================

nx, ny = size(∇f_x,2), size(∇f_y,2);
η_shocks  =  Ψ[ ny+1:end, : ]; ###  WARN:  select the equations refering to transition of states  ###
ϵ_sim = randn(1,100);

x_klein = zeros(_nVAR,101)
x_sims1 = zeros(_nVAR,101)
for it =1:100
    x_klein[1:nx,it+1] = hx * x_klein[1:nx,it] + η_shocks*ϵ_sim[:,it];
    x_klein[1+nx:end,it+1] = gx * x_klein[1:nx,it+1];
    x_sims1[:,it+1] = Phi1*x_sims1[:,it] + B1 * ϵ_sim[:,it];
    # x_sim2[:,it+1] = Phi2*x_sim2[:,it] + B2 * ϵ_sim[:,it];
end
maxabs(x_klein[:] - x_sims1[:])

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

# **************************************************************************************
#   FIGURES - simulation
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
using QuantEcon

H   = eye(size(Phi1,1));
lss = LSS(Phi1, B1, H)  ;# mu_0 = zeros(size(G, 2), Sigma_0 = zeros(size(H, 2), size(H, 2)) )

simul_length = 200
X_simul, _ = simulate(lss, simul_length);

Xss = [ss_histogram.xstst; ss_histogram.ystst];
X_simul = X_simul + repmat(Xss,1,simul_length);

nEps        = Cons.nEps;
nSavingsPar = Cons.nSavingsPar;
nAssetsFine = Cons.nAssetsFine;
asset_grid_fine = Cons.asset_grid_fine;
nHistogram  = nAssetsFine * Cons.nEps ;

K_simul            = X_simul[ss_histogram.ixvar[:aggr_end],: ];
z_simul            = X_simul[ss_histogram.ixvar[:aggr_exo],: ];

Θ_simul            = X_simul[ nx+ss_histogram.iyvar[:household],: ];
Cpol_simul         = zeros(nAssetsFine*nEps, simul_length);
xhistogram_simul      = X_simul[ ss_histogram.ixvar[:histogram],: ];
Histogram_simul = zeros(nAssetsFine,simul_length);

vHistogram  = zeros(nHistogram) ;
vHistogramK = zeros(nAssetsFine);
Cpol        = zeros(nAssetsFine*nEps);

for t = 1:simul_length

    #== Histogram ==#
    Reiter.x2distr!(vHistogram, xhistogram_simul[:,t])
    Histogram_simul[:,t] = sum( reshape(vHistogram, nAssetsFine, 2) , 2)

    #== Consumption ==#
    R      = 1 + Reiter.netintr(K_simul[t], z_simul[t]);
    wage   = Reiter.wagefunc(K_simul[t], z_simul[t]);
    for ieps = 1:nEps
        ind = (ieps-1)*nSavingsPar+1 : ieps*nSavingsPar
        Sthis= Reiter.interpolate(Cons, Θ_simul[ind,t])
        sthis = Sthis[asset_grid_fine]

        ind = (ieps-1)*nAssetsFine+1 : ieps*nAssetsFine
        Cpol[ind] = Reiter.get_consumption(asset_grid_fine, sthis, R, wage, ieps, Cons)
    end
    Cpol_simul[:,t] = Cpol
end

#== HIGH and LOW period ==#
ind = findmax(z_simul)[2], findmin(z_simul)[2];

Cpol_SS = zeros(nAssetsFine, nEps)
for ieps =1:nEps
    Sthis = Reiter.interpolate(Cons, ss_histogram.mΘ[:,ieps])
    sthis = Sthis[asset_grid_fine]
    Cpol_SS[:,ieps] = Reiter.get_consumption(asset_grid_fine, sthis, ss_histogram.R, ss_histogram.wage, ieps, Cons)
end

# .....................................................................................
using Plots, LaTeXStrings
plotlyjs()
y_vals = []
labels = String[]

push!(y_vals, Cpol_SS[:,1]); push!(labels, "stst unemp")
push!(y_vals, Cpol_simul[1:nAssetsFine, ind[1] ]); push!(labels, "High unemp ")
push!(y_vals, Cpol_simul[1:nAssetsFine, ind[2] ]); push!(labels, "Low unemp ")

push!(y_vals, Cpol_SS[:,2]); push!(labels, "Stst emp")
push!(y_vals, Cpol_simul[(1:nAssetsFine)+nAssetsFine, ind[1] ]); push!(labels, "High emp")
push!(y_vals, Cpol_simul[(1:nAssetsFine)+nAssetsFine, ind[2] ]); push!(labels, "Low  emp")

#== Consumption policies ==#
# plot(asset_grid_fine, y_vals[1], label=labels[1], linewidth=2, alpha=0.6, xlims=[0.0, 0.5*asset_grid_fine[end] ]);
psim1 = plot(asset_grid_fine, y_vals, label=labels', linewidth=2, alpha=0.6, xlims=(0.0, 6),xlabel="current assets", ylabel="Consumption", title =" Consumption policy");

#== Histogram ==#
palette = Plots.get_color_palette(:auto,default(:bgcolor),10)
psim2 = bar(Cons.asset_grid_fine, [Histogram_simul[:,ind[2]] Histogram_simul[:,1] Histogram_simul[:,ind[1]] ],
            label = ["low" "stst" "high"],
            color = [palette[2] palette[1] palette[3]],
            xlims=(0.0, 0.8*asset_grid_fine[end]), alpha=0.6,
            xlabel="Assets", title = "Distribution");

# COMMANDS to plot
plot(psim1, psim2, layout = 2, legendfont=font(10), guidefont = font(10))
plot(psim2)
# plot(p_hist_x_dens)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(-4, 4, 150)
y_vals = Array(Vector, 3)
labels = Array(String, 1, 3)
for i = 1:3
    m, s = 2*(rand() - 0.5), rand() + 1
    d = Normal(m, s)
    y_vals[i] = pdf(d, x)
    labels[i] = string(L"$\mu = $", round(m, 2))
end

plot(x, y_vals, linewidth=2, alpha=0.6, label=labels)

# ***************************************************************************************************************************************
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# ***************************************************************************************************************************************

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
