
using CompEcon
using NLsolve

fp = Reiter.FirmProblem(Spline,Cheb)

coeff = zeros(300,3);
ca = coeff[:,1];
cn = coeff[:,2];
cv = coeff[:,3];

#== Basis matrices ==#
Φ       = fp.mbasis.Φ;
Φ_z,_   = fp.mbasis.Φ_tensor.vals;

z_vals   = fp.mbasis.z_nodes;
grid_nodes = fp.mbasis.grid_nodes;

#== Indices ==#
ind_z_x_z = fp.mbasis.ind_z_x_z;
ind_z_x_p̃ = fp.mbasis.ind_z_x_p̃;

#==  ==#
p̃_basis = fp.mbasis.p̃_basis;

#== Guesses ==#
w, Y = 0.5, 1.0

###                               TEST Value Iteration                               ###
# **************************************************************************************
n_z     = fp.mbasis.n_z
n_p̃     = fp.mbasis.n_p̃
n_nodes = n_z*n_p̃

Φ_fac   = fp.mbasis.Φ_fac;
coeff = copy(fp.mbasis.coeff);
# .....................................................................................
tol = 1e-5
V    = fp.mbasis.Φ * coeff;
Vnew = similar(V);
ξstar= similar(V[:,1]);

ii = 1
err = 1.0
while ii <= 10 && err>tol
    Reiter.bellman_rhs!(Vnew, ξstar, coeff, w, Y, fp)

    # [Vnew ξstar]
    copy!(coeff,Φ_fac\Vnew)

    err = maximum( abs( vec(V) - vec(Vnew) ) )
    @printf("  value function iteration %d, distance %.6f \n", ii, err)

    #== Prepare NEXT iteration ==#
    copy!(V, Vnew)
    fill!(Vnew, 0.0)
    ii += 1
end

copy!(fp.coeff, coeff)
copy!(fp.ξstar, ξstar)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cv = coeff[:,3]
pstar = log( fp.ϵ/(fp.ϵ-1.0) * w./z_vals);

resid = zeros(10);
Reiter.foc_price_adjust!(resid, zeros(10), z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
resid

f!(p̃, fvec) = Reiter.foc_price_adjust!(fvec, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
g!(p̃, J)    = Reiter.soc_price_adjust!(J, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)

@time res1 = nlsolve(f!, zeros(10))
@time res2 = nlsolve(f!, g!, zeros(10))
[exp(res1.zero) exp(res2.zero) exp(pstar)]

#== FIND optimal p̃ in case of adjustment ==#

# **************************************************************************************
#   FIND optimal p̃ in case of adjustment
#
# ======================================================================================
# Solve for optimal price
pstar = log( fp.ϵ/(fp.ϵ-1.0) * w./z_vals);

resid = zeros(10);
Reiter.foc_price_adjust!(resid, zeros(10), z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
resid

f!(p̃, fvec) = Reiter.foc_price_adjust!(fvec, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
@time res = nlsolve(f!, zeros(10))
[exp(res.zero) exp(pstar)]
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# test function
resid = zeros(10);
Reiter.foc_price_adjust!(resid, zeros(10), z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
resid

resid[:] = 0.0;
Reiter.foc_price_adjust!(resid, pstar, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
resid

Jout = zeros(10,10);
Reiter.soc_price_adjust!(Jout, pstar, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
Jout

f!(p̃, fvec) = Reiter.foc_price_adjust!(fvec, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
g!(p̃, J)    = Reiter.soc_price_adjust!(J, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)

@time res = nlsolve(f!, zeros(10), autodiff = true)
@time res1 = nlsolve(f!, zeros(10))
@time res2 = nlsolve(f!, g!, zeros(10))
[exp(res1.zero) exp(res2.zero) exp(pstar)]
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
