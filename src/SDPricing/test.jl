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
w, Y = 1.0, 2.5

#== FIND optimal p̃ in case of adjustment ==#

# test function
resid = zeros(10)
Reiter.foc_price_adjust(resid, zeros(10), z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
resid


pstar = log( fp.ϵ/(fp.ϵ-1.0) * w./z_vals)
Reiter.foc_price_adjust(resid, pstar, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
resid

# Solve for optimal price
f!(p̃, fvec) = Reiter.foc_price_adjust(fvec, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)
g!(p̃, J)    = Reiter.soc_price_adjust(J, p̃, z_vals, w, Y, cv, fp, p̃_basis, Φ_z, ind_z_x_z)

@time res1 = nlsolve(f!, zeros(10))
@time res2 = nlsolve(f!,g!, ones(10))

exp([res1.zero res2.zero pstar])
