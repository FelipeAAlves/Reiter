


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


# x_t : k_t, z_t, y_t, c_t, i_t, n_t,
# x_ss = [ 3.49, 1.0, 0.90, 0.70, 0.21, 0.47]
np = -0.47/(1-0.47)
Γ0 = [
    -0.33 -1   0. 1 0 (0.33-1*np)       ;
    0.06 -0.08 0 1.5  0 (-0.06 + 0.5*np);
    -0.33 -1   1 0 0 -0.67              ;
    0 0        1 -0.77 -0.23 0          ;
    1 0        0 0 0 0                  ;
    0 1        0 0 0 0                  ]

Γ1 = [
    0 0            0. 0 0 0        ;
    0 0            0 1.5 0 (0.5*np);
    0 0            0 0 0 0         ;
    0 0            0 0 0 0         ;
    0.94 0         0 0 0.06 0      ;
    0 0.9          0 0 0 0         ]

Ψ = [0 0 0 0 0 1.0]'

Π = [0 1.0 0 0 0 0]'
c = zeros(6,1)



# using Gallium
using Gensys
include("C:\\Users\\falves\\Git\\Reiter\\src\\tools\\klein.jl")
# Gallium.breakpoint(gensysdt)

Phi, cons, B, _, _, _, _, eu, _ = gensysdt(Γ0, Γ1, c , Ψ , Π, 1.0 + 1e-10);
eu

Γ0 = [
    0 0            0. 0 0 0             ;
    0.06 -0.08 0 1.5  0 (-0.06 + 0.5*np);
    0 0             0 0 0 0             ;
    0 0             0 0 0 0             ;
    1 0        0 0 0 0                  ;
    0 1        0 0 0 0                  ]

Γ1 = [
    -0.33 -1   0. 1 0 (0.33-1*np)       ;
    0 0            0 1.5 0 (0.5*np);
    -0.33 -1   1 0 0 -0.67              ;
    0 0        1 -0.77 -0.23 0          ;
    0.94 0          0 0 0.06 0      ;
    0 0.9           0 0 0 0         ]

gx, hx = klein(Γ0[:,1:2], Γ0[:,3:end], -Γ1[:,1:2], -Γ1[:,3:end])
Phi, cons, B, _, _, _, _, eu, _ = gensysdt(Γ0, Γ1, c , Ψ , Π, 1.0 + 1e-10);
eu


nx,ny = 2,4
_nVAR = nx + ny
ϵ_sim = randn(10)

x_sim1 = zeros(_nVAR,11)
x_sim2 = zeros(_nVAR,11)
for it =1:10
    x_sim1[1:nx,it+1] = hx * x_sim1[1:nx,it] + [0;1]*ϵ_sim[it];
    x_sim1[nx+1:end,it+1] = gx * x_sim1[1:nx,it+1];
    x_sim2[:,it+1] = Phi*x_sim2[:,it] + B * ϵ_sim[it];
end
maxabs(x_sim1[:] - x_sim2[:])


    div = 1.0
    ϵ = 1e-6
    F = schurfact!(complex(Γ0), complex(Γ1))
    eu = [0, 0]
    a, b = F[:S], F[:T]

    n = size(a, 1)
    neta = size(Π, 2)

    for i in 1:n
        if (abs(a[i, i]) < ϵ) && (abs(b[i, i]) < ϵ)
            info("Coincident zeros.  Indeterminacy and/or nonexistence.")
            eu = [-2, -2]
            G1 = Array{Float64, 2}() ;  C = Array{Float64, 1}() ; impact = Array{Float64, 2}() ; fmat = Array{Complex{Float64}, 2}() ; fwt = Array{Complex{Float64}, 2}() ; ywt = Vector{Complex{Float64}}() ; gev = Vector{Complex{Float64}}() ; loose = Array{Float64, 2}()
            return G1, C, impact, fmat, fwt, ywt, gev, eu, loose
        end
    end

    movelast = Bool[abs(b[i, i]) > div * abs(a[i, i]) for i in 1:n]
    nunstab = sum(movelast)
    FS = ordschur!(F, !movelast)
    a, b, qt, z = FS[:S], FS[:T], FS[:Q], FS[:Z]


    gev = hcat(diag(a), diag(b))
    qt1 = qt[:, 1:(n - nunstab)]
    qt2 = qt[:, (n - nunstab + 1):n]
    q2xΠ = Ac_mul_B(qt2, Π)
    q2xΨ = Ac_mul_B(qt2, Ψ)
    q1xΠ = Ac_mul_B(qt1, Π)

    ndeta1 = min(n - nunstab, neta)

    # branch below is to handle case of no stable roots, rather than quitting with an error
    # in that case.
    if nunstab == 0
        q2xΠ = zeros(0, neta)
        ueta = zeros(0, 0)
        deta = zeros(0, 0)
        veta = zeros(neta, 0)
        bigev = 0
    else
        bigev, ueta, deta, veta = Gensys.decomposition_svd(q2xΠ)
        _    , teta, seta, weta = Gensys.decomposition_svd(q2xΨ)
    end

    # branch below to handle case of no stable roots
    if nunstab == n
        q1xΠ = zeros(0, neta)
        bigev1 = 0
        ueta1 = zeros(0, 0)
        veta1 = zeros(neta, 0)
        deta1 = zeros(0, 0)
    else
        bigev1, ueta1, deta1, veta1 = Gensys.decomposition_svd(q1xΠ)
    end


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# MATLAB

np = -0.47/(1-0.47)
Gamma0 = [
    -0.33 -1   0. 1 0 (0.33-1*np)       ;
    0.06 -0.08 0 1.5  0 (-0.06 + 0.5*np);
    -0.33 -1   1 0 0 -0.67              ;
    0 0        1 -0.77 -0.23 0          ;
    1 0        0 0 0 0                  ;
    0 1        0 0 0 0                  ];

Gamma1 = [
    0 0            0. 0 0 0        ;
    0 0            0 1.5 0 (0.5*np);
    0 0            0 0 0 0         ;
    0 0            0 0 0 0         ;
    0.94 0         0 0 0.06 0      ;
    0 0.9          0 0 0 0         ];

Psi = [0 0 0 0 0 1.0]';

Pi = [0 1.0 0 0 0 0]';
c = zeros(6,1);

[G1, C , impact, fmat, fwt, ywt, gev, eu1, loose]=gensys(Gamma0,Gamma1,c,Psi,Pi);


Gamma0 = [
    0 0             0. 0 0 0             ;
    0.06 -0.08      0 1.5  0 (-0.06 + 0.5*np);
    0 0             0 0 0 0             ;
    0 0             0 0 0 0             ;
    1 0             0 0 0 0                  ;
    0 1             0 0 0 0                  ]

Gamma1 = [
    -0.33 -1   0. 1 0 (0.33-1*np)       ;
    0 0            0 1.5 0 (0.5*np);
    -0.33 -1   1 0 0 -0.67              ;
    0 0        1 -0.77 -0.23 0          ;
    0.94 0          0 0 0.06 0      ;
    0 0.9           0 0 0 0         ];

[G1, C , impact, fmat, fwt, ywt, gev, eu2, loose]=gensys(Gamma0,Gamma1,c,Psi,Pi);

[eu1 eu2]
