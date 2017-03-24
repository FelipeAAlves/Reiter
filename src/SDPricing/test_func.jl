
# **************************************************************************************
#   USING code
#
# ======================================================================================
using BasisMatrices
using NLsolve
fcoll = Reiter.FirmColloc(Spline)
sol   = Reiter.FirmSolution(fcoll)

Reiter.solve_firm_policy!(fcoll, sol, 0.5, 1.0)

#== Finding the Steady-State ==#
ss_histogram = Reiter.StstHistogram(fcoll)

##  WARN : check result is not changing  ##
Reiter.stst_histogram_resid(0.75, 1.00, ss_histogram, fcoll, sol)
Reiter.stst_histogram_resid(1.25, 1.00,ss_histogram, fcoll, sol)


function f_hist!(x,fvec)
    copy!( fvec, Reiter.stst_histogram_resid(x[1], x[2], ss_histogram, fcoll, sol) )
end
res = nlsolve(f_hist!, [0.75,1.0/3], iterations = 50)

#== RUN one more time to UPDATE the ss_histogram ==#
Reiter.stst_histogram_resid(res.zero[1], res.zero[2], ss_histogram, fcoll, sol, true)

# **************************************************************************************
#   FIGURE : Histogram
#
# ======================================================================================

hist_nodes, (p_hist_nodes, z_hist_nodes) = Reiter.nodes(ss_histogram)
n_hist_p, n_hist_z = length(p_hist_nodes), length(z_hist_nodes)
vHistogram = Reiter.x2distr( ss_histogram.xstst[ss_histogram.ixvar[:histogram]] );
mHistogram = reshape(vHistogram, n_hist_p, n_hist_z) # p × z

#== pHistogram ==#
pHistogram = squeeze(sum(mHistogram,2),2)
pHistogram_plot = [sum(pHistogram[3*(i-1)+(1:3)]) for i=1:32]
Δp = p_hist_nodes[2]-p_hist_nodes[1];

#== In terms of deviation to optimal price ==#
pdev = exp( repmat(p_hist_nodes, n_hist_z) - kron(sol.pstar, ones(n_hist_p))   );
ind_dev = sortperm(pdev)
pdev = pdev[ind_dev]
vHistogram_dev = vHistogram[ind_dev]


p_x = linspace(pdev[1,1],pdev[end,1],101)
pdev_mass_plot = zeros(100)
i = 1
for (j,pp) in enumerate(p_x[2:end])
    while i<=length(pdev) && pdev[i,1]<=pp
        pdev_mass_plot[j] += vHistogram_dev[i]
        i+=1
    end
end
p_x = 0.5*(p_x[2:end]+p_x[1:end-1])

#== Δprice ==#
Φ_z   = fcoll.Φ_tensor.vals[2]
p_basis = fcoll.p_basis
cv = sol.coeff[:,3]
function eval_v̂(x::Vector{Float64}, z_ind::Vector{Int64}, deriv::Int64 = 0)

    @assert length(x)==length(z_ind)
    Φ_p̃ = BasisMatrix( p_basis, Direct(), x, deriv).vals[1]
    Φ_eval   = row_kron( Φ_z[ z_ind, :],  Φ_p̃ )
    return Φ_eval * cv
end

# get_xi_at{T<:Real}(p_nodes::Vector{Float64}, eval_v̂::Function, pᵃ::Array{Float64}, w::T, Y::T, Y_pri::T=Y, Π_pri::T = 1.0, z_aggr::T=0.0)
ξstar_distr, _ = Reiter.get_xi_at(p_hist_nodes, eval_v̂, sol.pstar, ss_histogram.ystst[5], ss_histogram.ystst[1]);
Δp_distr_mass = vHistogram .* Reiter.__pars.H(ξstar_distr);
Δp_distr = exp( kron(sol.pstar, ones(n_hist_p)) - repmat(p_hist_nodes, n_hist_z)  );

ind_Δ = sortperm(Δp_distr)
Δp = [Δp_distr Δp_distr_mass][ind_Δ,:]
Δp_x = linspace(Δp[1,1],Δp[end,1],81)
Δp_mass_plot = zeros(80)

i = 1
for (j,pp) in enumerate(Δp_x[2:end])
    while i<=length(Δp_distr) && Δp[i,1]<=pp
        Δp_mass_plot[j] += Δp[i,2]
        i+=1
    end
end
Δp_x = 0.5*(Δp_x[2:end]+Δp_x[1:end-1])

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using Plots
# using LaTeXStrings
# pyplot() #
plotlyjs()

#== PLOT 01 ==#
title = "Stationary Distribution of Prices"
plt_stat_dist =  bar(exp(p_hist_nodes),  pHistogram ./ (Δp), alpha=0.6, label="psi",  legendfont=font(12));
plot(plt_stat_dist, title=title) # in order to show

#== PLOT 02 ==#
plot( Vector{Float64}[, ], Vector{Float64}[Δp_mass_plot./sum( Δp_mass_plot ),  pHistogram], alpha=0.6, label="psi",  legendfont=font(12))
bar( Vector{Float64}[1:10, 1:0.5:10], Vector{Float64}[rand(10),  rand(19)], alpha=0.6, label="psi",  legendfont=font(12), Layout(barmode="overlay"))
bar( Vector{Float64}[0.5*(Δp_x[2:end]+Δp_x[1:end-1]), exp(p_hist_nodes)], Vector{Float64}[Δp_mass_plot./sum( Δp_mass_plot ),  pHistogram], alpha=0.6, label="psi",  legendfont=font(12))
plot!(exp(p_hist_nodes),  pHistogram, alpha=0.6, label="psi",  legendfont=font(12))
plt_Δp = bar(0.5* (Δp_x[2:end]+Δp_x[1:end-1]), Δp_mass_plot ./ sum( Δp_mass_plot ), alpha=0.6, label="psi",  legendfont=font(12))
plot(plt_stat_dist,plt_Δp)
plot!(plt_Δp,xaxis=("Price Change",(0.55,1.45),0.6:0.1:1.40))

# .....................................................................................
# .....................................................................................

using PlotlyJS
x = linspace(0, 10, 200)
y = sin(x)
plot(scatter(x=x, y=y, marker_color="blue", line_width=2))

function two_distributions()
    trace1 = bar(;x=Δp_x,y=Δp_mass_plot./sum( Δp_mass_plot ),opacity=0.75, name="Price Changes")
    trace2 = bar(;x=p_x, y=pdev_mass_plot,opacity=0.75, name="Price Changes")
    # trace2 = bar(;x=exp(p_hist_nodes[2:3:end-1]),y=pHistogram_plot,opacity=0.75, name="Stationary Distribution")
    # trace2 = bar(;x=exp(p_hist_nodes),y=pHistogram,opacity=0.75, name="Stationary Distribution")

    data = [trace1, trace2]
    #== Separate ==#
    layout1 = Layout(;barmode="overlay",title="Stationary Distribution",xaxis_range=[0.55,1.45], xaxis_title="Deviation"); p1 = plot(data[1], layout1)
    layout2 = Layout(;barmode="overlay",title="Price Changes",xaxis_range=[0.55,1.45], xaxis_title="Deviation"); p2 = plot(data[2], layout2)
    [p1 p2]

    #== TOGETHER ==#
    # layout = Layout(;barmode="overlay",title="Stationary Distribution",xaxis_range=[0.55,1.45], xaxis_title="Deviation");
    # plot(data, layout1)
end
two_distributions()


# **************************************************************************************
#   LINEARIZING
#
# ======================================================================================

res_equil = Reiter.equil_histogram(ss_histogram.Zstst, ss_histogram, fcoll, sol);
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

f(Y) = Reiter.equil_histogram(Y, ss_histogram, fcoll, sol)
∇f = zeros( length(res_equil), length(ss_histogram.Zstst) );
@time ForwardDiff.jacobian!(∇f, f, ss_histogram.Zstst, ForwardDiff.JacobianConfig{5}(ss_histogram.Zstst));
# @time ∇f = ForwardDiff.jacobian(f, Zss, ForwardDiff.JacobianConfig{5}(Zss))


#================================================================#
###                  Applying Klein algorithm                  ###
#================================================================#
∇f_y′   = ∇f[:,  ss_histogram.iZvar[:y′] ] ;
∇f_x′   = ∇f[:,  ss_histogram.iZvar[:x′] ] ;
∇f_y    =  ∇f[ :, ss_histogram.iZvar[:y] ]  ;
∇f_x    =  ∇f[ :, ss_histogram.iZvar[:x] ]  ;
η_shocks    =  -∇f[ length(ss_histogram.ystst)+1:end, ss_histogram.iZvar[:eps] ]; # equations refering to transition of states

gx, hx = Reiter.klein(∇f_x′, ∇f_y′, ∇f_x, ∇f_y);

# **************************************************************************************
# **************************************************************************************

Γ0 = -[ ∇f_x′ ∇f_y′ ];
Γ1 =  [ ∇f_x ∇f_y ];
c  = zeros(_nVAR);
Ψ  = ∇f[ :, ss_histogram.iZvar[:eps] ];

_Pi_  = - ∇f_y′; # E_t{ y_{t+1} } = y_t+1 - eta_t+1
ind_η = Reiter.findnzcols(_Pi_, 1e-12) # Reiter.findnzcols(_Pi_);
_nPi  = length(ind_η)
Π     = _Pi_[:, ind_η ]

#== effective/contaminated foward looking equations ==#
# ind_Pi = Reiter.findnzrows(_Pi_)
# ind_Pi_effec = ind_Pi[1:_nPi]
# ind_Pi_extra = ind_Pi[_nPi+1:end]
#
# Λ = _Pi_[ind_Pi_effec,:]
# _Pi = spzeros(_nVAR,_nPi)
# _Pi[ind_Pi_effec,:] = speye(_nPi)
# _Pi[ind_Pi_extra,:] = _Pi_[ind_Pi_extra,ind_η]/Λ

# using Gensys
@time Phi, cons, B, _, _, _, _, eu, _ = Reiter.gensysdt(Γ0, Γ1, c , Ψ , Π, 1.0 + 1e-10);
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



# ***************************************************************************************************************************************
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# ***************************************************************************************************************************************


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


#===========================================================#
##                        MAIN usage                       ##
#===========================================================#
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


# reshape(ss_histogram.Xss[ss_histogram.iXvar[:value_fnc]],32,10)

# .....................................................................................

# using Gallium
# Gallium.breakpoint(Reiter.equil_histogram)
