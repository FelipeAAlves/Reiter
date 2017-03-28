
# **************************************************************************************
#   USING code
#
# ======================================================================================
using BasisMatrices
using NLsolve
fcoll = Reiter.FirmColloc(Spline); exp(fcoll.p_nodes)
sol   = Reiter.FirmSolution(fcoll)

Reiter.solve_firm_policy!(fcoll, sol, 0.80, 0.33)
sol.ξstar
#== Finding the Steady-State ==#
ss_histogram = Reiter.StstHistogram(fcoll)

##  WARN : check result is not changing  ##

Reiter.stst_histogram_resid(0.65, 0.33, ss_histogram, fcoll, sol)
Reiter.stst_histogram_resid(0.75, 0.33, ss_histogram, fcoll, sol)
Reiter.stst_histogram_resid(0.85, 0.33, ss_histogram, fcoll, sol)
Reiter.stst_histogram_resid(0.95, 0.33, ss_histogram, fcoll, sol)

Reiter.stst_histogram_resid(0.800, 0.32, ss_histogram, fcoll, sol)
function f_hist!(x,fvec)
    copy!( fvec, Reiter.stst_histogram_resid(x[1], x[2], ss_histogram, fcoll, sol) )
end
res = nlsolve(f_hist!, [0.800, 0.32], iterations = 50)

# res = nlsolve(f_hist!, [0.75, 0.32], iterations = 10, method = :newton) # NOTE worst option
# res = nlsolve(f_hist!, [0.80,Reiter.__pars.Nstst], iterations = 100, autodiff = true)

#== RUN one more time to UPDATE the ss_histogram ==#
Reiter.stst_histogram_resid(res.zero[1], res.zero[2], ss_histogram, fcoll, sol, true)

# **************************************************************************************
#   FIGURE : Histogram
#
# ======================================================================================
# if plot
    hist_nodes, (p_hist_nodes, z_hist_nodes) = Reiter.nodes(ss_histogram)
    n_hist_p, n_hist_z = length(p_hist_nodes), length(z_hist_nodes)

    #== Effective prices distribution ==#
    Φ_z  = fcoll.Φ_tensor.vals[2]
    p_basis = fcoll.p_basis
    cv = sol.coeff[:,3]
    function eval_v̂(x, z_ind::Vector{Int64}, deriv::Int64 = 0)

        @assert length(x)==length(z_ind)
        Φ_p̃ = BasisMatrix( p_basis, Direct(), x, deriv).vals[1]
        Φ_eval   = row_kron( Φ_z[ z_ind, :],  Φ_p̃ )
        return Φ_eval * cv
    end
    ξstar_distr, _ = Reiter.get_xi_at(p_hist_nodes, eval_v̂, sol.pstar, res.zero[1], res.zero[2])
    Π1 = Reiter.endperiod_transition( p_hist_nodes, sol.pstar, ξstar_distr; update_idio=false ) ###  WARN:  NO-adjustment of idio values  ###
    vHistogram = Reiter.x2distr( ss_histogram.xstst[ss_histogram.ixvar[:histogram]] );
    vHistogram_end_of_period = Π1 * vHistogram;
    mHistogram = reshape(vHistogram_end_of_period, n_hist_p, n_hist_z) # p × z

    #== pHistogram ==#
    pHistogram = squeeze(sum(mHistogram,2),2)
    pHistogram_plot = [sum(pHistogram[3*(i-1)+(1:3)]) for i=1:32]
    Δp = p_hist_nodes[2]-p_hist_nodes[1];

    #== In terms of relative price ==#
    phist_x = linspace(exp(p_hist_nodes[1]),exp(p_hist_nodes[end]),101)
    phist_mass_plot = zeros(100)
    i = 1
    for (j,pp) in enumerate(phist_x[2:end])
        while i<=length(p_hist_nodes) && exp(p_hist_nodes[i])<=pp
            phist_mass_plot[j] += pHistogram[i]
            i+=1
        end
    end
    phist_x = 0.5*(phist_x[2:end]+phist_x[1:end-1])


    #== In terms of deviation to optimal price ==#
    pdev = exp( repmat(p_hist_nodes, n_hist_z) - kron(sol.pstar, ones(n_hist_p))   );
    ind_dev = sortperm(pdev)
    pdev = pdev[ind_dev]
    vHistogram_dev = vHistogram_end_of_period[ind_dev]


    pdev_x = linspace(pdev[1],pdev[end],101)
    pdev_mass_plot = zeros(100)
    i = 1
    for (j,pp) in enumerate(pdev_x[2:end])
        while i<=length(pdev) && pdev[i]<=pp
            pdev_mass_plot[j] += vHistogram_dev[i]
            i+=1
        end
    end
    pdev_x = 0.5*(pdev_x[2:end]+pdev_x[1:end-1])

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

    ξstar_distr, _ = Reiter.get_xi_at(p_hist_nodes, eval_v̂, sol.pstar, ss_histogram.ystst[5], ss_histogram.ystst[1]);
    Δp_distr_mass = vHistogram_end_of_period .* Reiter.__pars.H(ξstar_distr);
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
    plotlyjs()
    #== TEST 01 ==#
    x = linspace(0, 10, 100)
    plot(x,sin,color=:red,lw=2,yticks=-1:1:1,title="sine function",alpha=0.6)

    #== TEST 02 ==#
    x = linspace(-4, 4, 150)
    y_vals = Array(Vector, 3)
    labels = Array(String, 1, 3)
    for i = 1:3
        m, s = 2*(rand() - 0.5), rand() + 1
        d = Normal(m, s)
        y_vals[i] = pdf(d, x)
        labels[i] = string("mu = ", round(m, 2))
    end

    plot(x, y_vals, linewidth=2, alpha=0.6, label=labels)
    plot(Plots.fakedata(100,10),layout=4,palette=[:grays :blues :heat :lightrainbow],bg_inside=[:orange :pink :darkblue :black])

    #== PLOT 01 ==#
    title = "Stationary Distribution of Prices"
    plt_stat_dist =  bar( phist_x,  phist_mass_plot, alpha=0.6, label="Stationary Distribution",  legendfont=font(12) );
    plot(plt_stat_dist, title=title) # in order to show

    #== PLOT 02 ==#
    # test
    palette = Plots.get_color_palette(:auto,default(:bgcolor),10)
    bar(Vector{Float64}[1:10, 1:0.5:10], Vector{Float64}[rand(10),  rand(19)], layout=2, alpha=0.6, color=[palette[1] palette[2]], label=["psi1","psi2"],  legendfont=font(10), xaxis=("Relative price",(1,12),1:0.5:10))
    bar(1:10, [rand(10)  rand(10)], alpha=0.6, label=["psi1","psi2"],  legendfont=font(10), xaxis=("Relative price",(1,12),1:0.5:10))

    # getcol(i) = Plots.getSeriesRGBColor(:auto,Dict(:color_palette=>palette),i)
    bar( Vector{Float64}[pdev_x, Δp_x], Vector{Float64}[pdev_mass_plot, Δp_mass_plot],
            alpha=0.6, label=["Relative prices (p-pstar)","Price Changes"],
            xaxis=("Relative price",(0.55,1.45),0.6:0.1:1.40), color = [getcol(1) getcol(2)],
            legendfont=font(10), layout =2 )
    plot!(plt_Δp,xaxis=("Price Change",(0.55,1.45),0.6:0.1:1.40))

    # .....................................................................................
    # .....................................................................................

    using PlotlyJS
    x = linspace(0, 10, 200)
    y = sin(x)
    plot(scatter(x=x, y=y, marker_color="blue", line_width=2))

    function two_distributions()
        trace0 = bar(;x=phist_x, y=phist_mass_plot,opacity=0.75, name="Relative prices (stationary dist)")
        trace1 = bar(;x=pdev_x, y=pdev_mass_plot,opacity=0.75, name="Relative prices (p-p_star)")
        trace2 = bar(;x=Δp_x,y=Δp_mass_plot./sum( Δp_mass_plot ),opacity=0.75, name="Price Changes")
        # trace2 = bar(;x=exp(p_hist_nodes[2:3:end-1]),y=pHistogram_plot,opacity=0.75, name="Stationary Distribution")
        # trace2 = bar(;x=exp(p_hist_nodes),y=pHistogram,opacity=0.75, name="Stationary Distribution")

        #== OPTION Separate ==#
        # layout1 = Layout(;barmode="overlay",title="Stationary Distribution",xaxis_range=[0.55,1.45], xaxis_title="Deviation"); p1 = plot(trace1, layout1)
        # layout2 = Layout(;barmode="overlay",title="Price Changes",xaxis_range=[0.55,1.45], xaxis_title="Deviation"); p2 = plot(trace2, layout2)
        # [p1 p2]

        #== OPTION same graph ==#
        data = [trace0, trace1, trace2]
        layout = Layout(;barmode="overlay",title="Stationary Distribution", xaxis_range=[0.55,1.45], xaxis_title="Deviation");
        plot(data, layout)
    end
    two_distributions()
end

# **************************************************************************************
#   LINEARIZING
#
# ======================================================================================
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

res_equil = Reiter.equil_histogram(ss_histogram.Zstst, ss_histogram, fcoll, sol);
maximum( abs( res_equil) ) > 1e-6  && throw(error("Wrong Stst"))
_nVAR = length(res_equil)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
using ForwardDiff

f(Y) = Reiter.equil_histogram(Y, ss_histogram, fcoll, sol)
∇f = zeros( length(res_equil), length(ss_histogram.Zstst) );
@time ForwardDiff.jacobian!(∇f, f, ss_histogram.Zstst, ForwardDiff.JacobianConfig{5}(ss_histogram.Zstst));
# @time ∇f = ForwardDiff.jacobian(f, Zss, ForwardDiff.JacobianConfig{5}(Zss))

∇f_x′   = ∇f[:,  ss_histogram.iZvar[:x′] ] ;
∇f_y′   = ∇f[:,  ss_histogram.iZvar[:y′] ] ;
∇f_x    =  ∇f[ :, ss_histogram.iZvar[:x] ]  ;
∇f_y    =  ∇f[ :, ss_histogram.iZvar[:y] ]  ;


#--------------------------------#
##             KLEIN             ##
#--------------------------------#
#== CHECK stable/unstable eigenvalues ==#
Γ0 =  [  ∇f_x′ ∇f_y′];
writecsv("gamma0.csv", Γ0)
Γ1 = -[ ∇f_x ∇f_y ];
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

Γ0 =  [ ∇f_x′ ∇f_y′];
Γ1 = -[ ∇f_x ∇f_y];
c  = zeros(_nVAR);
Ψ  = -∇f[ :, ss_histogram.iZvar[:eps] ];

_Pi_  = ∇f_y′; # E_t{ y_{t+1} } = y_t+1 - eta_t+1
    #== OPTION 01 ==#
    # works well... it is more natural way to do it I think
    # associate one expectational error to each control that appears on t+1
    # note, from the way I wrote the model, all controls appearing at t+1 are E_t\{\}

    ind_η = Reiter.findnzcols(_Pi_, 1e-12) # Reiter.findnzcols(_Pi_);
    _nPi  = length(ind_η)
    Π1     = _Pi_[:, ind_η ]

    #== OPTION 02 ==#
    ind_η2 = Reiter.findnzrows(_Pi_,1e-12)
    Π2     = zeros(_nVAR, length(ind_η2))
    for i=1:length(ind_η2)
        Π2[ind_η2[i],i]     = 1.0
    end

    writecsv("_Pi_.csv", _Pi_)
    writecsv("_Pi_1.csv", Π1)
    writecsv("_Pi_2.csv", Π2)

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
@time Phi1, cons, B1, _, _, _, _, eu1, _ = Reiter.gensysdt(Γ0, Γ1, c , Ψ , Π1, 1.0 + 1e-10);
@time Phi2, cons, B2, _, _, _, _, eu2, _ = Reiter.gensysdt(Γ0, Γ1, c , Ψ , Π2, 1.0 + 1e-10);
eu1
eu2

# **************************************************************************************
#   SIMULATION
#
# ======================================================================================
nx,ny = size(∇f_x,2), size(∇f_y,2)
η_shocks =  Ψ[ny+1:end, : ]; ###  WARN:  select equations referring to state transition     ###
ϵ_sim = randn(2,100)

x_klein = zeros(_nVAR,101)
x_sims1 = zeros(_nVAR,101)
x_sims2 = zeros(_nVAR,101)
for it =1:100
    x_klein[1:nx,it+1] = hx * x_klein[1:nx,it] + η_shocks*ϵ_sim[:,it];
    x_klein[1+nx:end,it+1] = gx * x_klein[1:nx,it+1];
    x_sims1[:,it+1] = Phi1*x_sims1[:,it] + B1 * ϵ_sim[:,it];
    x_sims2[:,it+1] = Phi2*x_sims2[:,it] + B2 * ϵ_sim[:,it];
end
maxabs(x_klein[:] - x_sims1[:])
maxabs(x_sims1[:] - x_sims2[:])



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
