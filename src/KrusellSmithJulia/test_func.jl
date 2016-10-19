

using Plots
pyplot()
using LaTeXStrings


using NLsolve
using Roots
using ForwardDiff

# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------
#== Initialize Objects ==#
Cons         = Reiter.ConsumerProblem();
ss_histogram = Reiter.StstHistogram(Cons);
#== Initialize Objects ==#
Cons         = Reiter.ConsumerProblem();
ss_histogram = Reiter.StstHistogram(Cons);

#== GUESS prices ==#
Rss      = 1 + Reiter.netintr(1.01*Cons.KrepSS, 0);
wagess   = Reiter.wagefunc(1.01*Cons.KrepSS, 0);

##  TEST eulerres fnc  ##
out1 = zeros(120*2)
out2 = zeros(120*2)
@time Reiter.eulerres1!(out1, Cons.Θinit , Cons.Θinit, Rss, Rss, wagess, wagess, Cons)
@time Reiter.eulerres2!(out2, Cons.Θinit , Cons.Θinit, Rss, Rss, wagess, wagess, Cons)

@code_warntype Reiter.eulerres1!(out1, Cons.Θinit , Cons.Θinit, Rss, Rss, wagess, wagess, Cons)
@code_warntype Reiter.eulerres2!(out2, Cons.Θinit , Cons.Θinit, Rss, Rss, wagess, wagess, Cons)

##  TEST policy function  ##
f!(Θ, fvec) = Reiter.eulerres1!(fvec, Θ , Θ, Rss, Rss, wagess, wagess, Cons)
g!(Θ, fvec) = Reiter.eulerres2!(fvec, Θ , Θ, Rss, Rss, wagess, wagess, Cons)

#== Compare performance ==#
# @time nlsolve(f!, Cons.Θinit)
@time nlsolve(f!, Cons.Θinit, autodiff = true)
@time nlsolve(g!, Cons.Θinit, autodiff = true)

res1 = nlsolve(f!, Cons.Θinit, autodiff = true)
res2 = nlsolve(g!, Cons.Θinit, autodiff = true)
copy!(Cons.Θinit, res2.zero);
# ..............................................................................

# ******************************************************************
#   TEST Derivatives
# ==================================================================
using ForwardDiff

## TEST 01 : test the jacobian on the INTEPOLATION ##
out = zeros(120*2);
f!(y, Θ) = Reiter.test_jaco_interp!(y, Θ, Rss, wagess, Cons)
J = ForwardDiff.jacobian( f!, out, Cons.Θinit, Chunk{10}())

@time ForwardDiff.jacobian( f!, out, Cons.Θinit, Chunk{2}() );
@time ForwardDiff.jacobian( f!, out, Cons.Θinit, Chunk{5}() );
@time ForwardDiff.jacobian( f!, out, Cons.Θinit, Chunk{10}());

##  TEST 02 : Do DUAL by hand  ##
Θ     = Dual{2,Float64}[Dual(Cons.Θinit[i],i==1,i==2 ) for i=1:240]
Θnext = Dual{2,Float64}[Dual(Cons.Θinit[i],0,0) for i=1:240]

num, cons = Reiter.eulerres_test(Θ, Θnext, Rss, Rss, wagess, wagess, Cons)
# CHECK the value with what comes out of jacobian function
1-num[2]/cons[2]

f!(y, Θall) = Reiter.eulerres2!(y, Θall[1:240], Θall[241:end], Rss, Rss, wagess, wagess, Cons)
@time J = ForwardDiff.jacobian(f!, out, repmat(Cons.Θinit,2) );
# .....................................................................................

## TEST 03 : test the jacobian on Euler Residuals ##
f!(y, Θall) = Reiter.eulerres2!(y, Θall[1:240], Θall[241:end], Rss, Rss, wagess, wagess, Cons)

#== No pre-allocation ==#
@time J = ForwardDiff.jacobian(f!, out, repmat(Cons.Θinit,2) );

#== Pre-allocated ==#
Jout = ForwardDiff.JacobianResult(out);
@time ForwardDiff.jacobian!(Jout, f!, Cons.Θinit, Chunk{10}());
@time ForwardDiff.jacobian!(Jout, f!, Cons.Θinit, Chunk{5}());
@time ForwardDiff.jacobian!(Jout, f!, Cons.Θinit, Chunk{2}());

Jout = similar(J)
f(Θ) = eulerres(Cons.Θinit, Θ, Rss, Rss, wagess, wagess, Cons)
@time J = ForwardDiff.jacobian!(Jout, f, Cons.Θinit);
# .....................................................................................

# ******************************************************************
#   TEST Histogram SteadyState
# ==================================================================

#== Initialize Objects ==#
Cons         = Reiter.ConsumerProblem();
ss_histogram = Reiter.StstHistogram(Cons);

f(k) = Reiter.stst_histogram_resid(Cons, k)
Kss = fzero(f, 1.01*Cons.KrepSS, 1.05*Cons.KrepSS)

#== UPDATE ss_histogram ==#
Reiter.stst_histogram_resid(Cons, Kss, ss_histogram)

# ******************************************************************
#   TEST Density SteadyState
# ==================================================================

ss_density = Reiter.StstDensity(Cons);

#== Compute moments from histogram  ==#
mMoments = zeros(Cons.nMoments, Cons.nEps);
for ieps = 1:Cons.nEps
    for iMom = 1:Cons.nMoments
        dens = ss_histogram.mHistogram[:,ieps]/sum( ss_histogram.mHistogram[:,ieps] )
        mMoments[iMom, ieps] = dot( (Cons.asset_grid_fine - (iMom>1) * mMoments[1, ieps] ).^iMom, dens )
    end
end

f(k) = Reiter.stst_density_resid(Cons, k, mMoments)
Kss_density = fzero(f, 0.995*ss_histogram.Kaggr, 1.005*ss_histogram.Kaggr)

#== UPDATE ss_histogram ==#
Reiter.stst_density_resid(Cons, Kss_density, mMoments, ss_density)

# ******************************************************************
#   FIGURE : Histogram × Density
# ==================================================================
using Plots
pyplot() # plotlyjs()
using LaTeXStrings

asset_grid_fine = Cons.asset_grid_fine;
asset_quad_points = Cons.asset_quad_points;
Δa = asset_grid_fine[2]-asset_grid_fine[1];
z_dist = Cons.z_dist

#== CHECK density ==#
Reiter.density_int(ss_density.mDensityCoeff[2:end,1], ss_density.mMoments[:,1], Cons, 0, ss_density.mDensityCoeff[1,1])
Reiter.density_int(ss_density.mDensityCoeff[2:end,2], ss_density.mMoments[:,2], Cons, 0, ss_density.mDensityCoeff[1,2])

ieps = 1
    p1 = bar(asset_grid_fine, ss_histogram.mHistogram[:,ieps]/(Δa*z_dist[ieps]), alpha=0.6, label=LaTeXString("\$ h_1(a_i)\$"),  legendfont=font(12));
    dens = ss_density.mDensityCoeff[1,ieps] * Reiter.density(ss_density.mDensityCoeff[2:end,ieps], ss_density.mMoments[:,ieps], asset_grid_fine, Cons);
    plot!(asset_grid_fine, dens, color=:red, linewidth=2, label=LaTeXString("\$ g_{1} (a) \$"),  legendfont=font(12));
ieps = 2
    p2 = bar(asset_grid_fine, ss_histogram.mHistogram[:,ieps]/(Δa*z_dist[ieps]), alpha=0.6,label=LaTeXString("\$ h_2(a_i, 1_j)\$"),  legendfont=font(12));
    dens = ss_density.mDensityCoeff[1,ieps] * Reiter.density(ss_density.mDensityCoeff[2:end,ieps], ss_density.mMoments[:,ieps], asset_grid_fine, Cons);
    plot!(asset_grid_fine, dens, color=:red, linewidth=2, label=LaTeXString("\$g_{2}(a)\$"), legendfont=font(12));

titles = ["Unemp" "Employed"]
plot(p1, p2, layout=2, title=titles)

# ******************************************************************
#   LINEARIZATION Step
# ==================================================================
println("STEP [2]: Linearizing model equations")
println("STEP [3]: Solving the model")

###  Construct Y  ###
# *************************************************************************************
#== Construct iYvar==#
iYvar = Dict{Symbol, Union{UnitRange{Int64}, Int64}}()

nEta = Cons.nSavingsPar * Cons.nEps
nExog = 1

#== CREATE Y ==#
# X = [vHistogramDev; Kaggr; dx; vSavingsPar]
Xss = [ zeros(length(ss_histogram.vHistogram)-1); 0.0; ss_histogram.Kaggr; ss_histogram.mΘ[:] ];
Yss = [ Xss; Xss; zeros(nEta); zeros(nExog) ];

#== iYvar ==#
nX = length(Xss)
iYvar[:X]    = 1:nX
iYvar[:Xlag] = nX + (iYvar[:X])

nn = 2*nX
iYvar[:eta] = nn+1:nn+nEta; nn = nn+nEta
iYvar[:eps] = nn+1
# .....................................................................................

#== Evaluate Stst ==#
equil_res = Reiter.equil_histogram(Yss, iYvar, ss_histogram, Cons)
maximum( abs( equil_res) ) > 1e-6  && throw(error("Wrong Stst"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
using ForwardDiff
using Gensys
using QuantEcon

f(Y) = Reiter.equil_histogram(Y, iYvar, ss_histogram, Cons)
Jout = zeros( length(equil_res), length(Yss) );
ForwardDiff.jacobian!(Jout, f, Yss);
# @time ForwardDiff.jacobian!(Jout, f, Yss,Chunk{10}());
# @time ForwardDiff.jacobian!(Jout, f, Yss,Chunk{5}()) ;
# @time ForwardDiff.jacobian!(Jout, f, Yss,Chunk{2}()) ;

#== Gensys matrices ==#
Γ0 = -Jout[:, iYvar[:X]]    ;
Γ1 =  Jout[:, iYvar[:Xlag]] ;
c  =  zeros( size(Jout,1) ) ;
Π  =  Jout[:, iYvar[:eta]]  ;
Ψ  =  Jout[:, iYvar[:eps]]; Ψ = reshape(Ψ, length(Ψ), 1 ) ;

#== CALL Gensys ==#
# ```
# Γ0*y(t) = Γ1*y(t-1) + c + Ψ*z(t) + Π*η(t),
# ```
Phi, cons, B, _, _, _, _, eu, _ = gensysdt(Γ0, Γ1, c, Ψ, Π, 1.0 + 1e-10);
if !((eu[1] == 1) & (eu[2] == 1))
    throw(error("Gensys does not give existence"))
end

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


# ******************************************************************
#   FIGURES
# ==================================================================
H   = eye(size(Phi,1));
lss = LSS(Phi, B, H)  ;# mu_0 = zeros(size(G, 2), Sigma_0 = zeros(size(H, 2), size(H, 2)) )

simul_length = 200
X_simul, _ = simulate(lss, simul_length);
X_simul = X_simul + repmat(Xss,1,simul_length);

nEps        = Cons.nEps
nSavingsPar = Cons.nSavingsPar
nAssetsFine = Cons.nAssetsFine
nHistogram  = nAssetsFine * Cons.nEps
iXvar              = ss_histogram.iXvar

K_simul            = X_simul[ iXvar[:aggr_end],: ]';
z_simul            = X_simul[ iXvar[:aggr_exo],: ]';

Θ_simul            = X_simul[ iXvar[:household],: ];
Cpol_simul         = zeros(nAssetsFine*nEps, simul_length);
HistogramDev_simul = X_simul[ iXvar[:histogram],: ];
Histogram_simul    = zeros(nAssetsFine,simul_length);

vHistogram  = zeros(nHistogram) ;
vHistogramK = zeros(nAssetsFine);
Cpol_simul  = zeros(nAssetsFine*nEps, simul_length);

asset_grid_fine = Cons.asset_grid_fine
for t = 1:simul_length

    #== Histogram ==#
    Reiter.x2distr!(vHistogram, HistogramDev_simul[:,t], ss_histogram.vHistogram)
    copy!(vHistogramK, sum( reshape(vHistogram, nAssetsFine, 2) , 2) )
    Histogram_simul[:,t] = vHistogramK

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

#== HIGH and Low period ==#
ind = findmax(z_simul)[2], findmin(z_simul)[2];

Cpol_SS = zeros(nAssetsFine, nEps)
for ieps =1:nEps
    Sthis = Reiter.interpolate(Cons, ss_histogram.mΘ[:,ieps])
    sthis = Sthis[asset_grid_fine]
    Cpol_SS[:,ieps] = Reiter.get_consumption(asset_grid_fine, sthis, ss_histogram.R, ss_histogram.wage, ieps, Cons)
end

# .....................................................................................

using Plots
using LaTeXStrings        # Install this package
# plotlyjs()
pyplot()

y_vals = []
labels = []

push!(y_vals, Cpol_SS[:,1])
push!(y_vals, Cpol_simul[1:nAssetsFine, ind[1] ])
push!(y_vals, Cpol_simul[1:nAssetsFine, ind[2] ])

push!(y_vals, Cpol_SS[:,2])
push!(y_vals, Cpol_simul[(1:nAssetsFine)+nAssetsFine, ind[1] ])
push!(y_vals, Cpol_simul[(1:nAssetsFine)+nAssetsFine, ind[2] ])

plot(asset_grid_fine, y_vals, linewidth=2, alpha=0.6, xlim=[0.0, 0.5*asset_grid_fine[end] ])

plot(asset_grid_fine, Cpol_SS[:,1], linewidth=2, alpha=0.6)

#== Histogram ==#
bar(Cons.asset_grid_fine, [Histogram_simul[:,1] Histogram_simul[:,ind[1]] Histogram_simul[:,ind[2]]], alpha=0.6)
# -------------------------------------------------------------------------------------