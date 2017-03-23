# This code is based on a routine originally copyright Chris Sims.
# See http://sims.princeton.edu/yftp/gensysdt/


"""
Options to call
```_
gensysdt(Γ0, Γ1, c, Ψ, Π)
gensysdt(Γ0, Γ1, c, Ψ, Π, div)
gensysdt(F::Base.LinAlg.GeneralizedSchur, c, Ψ, Π)
gensysdt(F::Base.LinAlg.GeneralizedSchur, c, Ψ, Π, div)
```

Generate state-space solution to canonical-form DSGE model.

System given as
```_
Γ0*y(t) = Γ1*y(t-1) + c + Ψ*z(t) + Π*η(t),
```
with z an exogenous variable process and η being endogenously
determined one-step-ahead expectational errors.

Returned system is
```_
y(t) = G1*y(t-1) + C + impact*z(t) + ywt*inv(I-fmat*inv(L))*fwt*z(t+1)
```
Returned values are
```_
G1, C, impact, fmat, fwt, ywt, gev, eu, loose
```

If `z(t)` is i.i.d., the last term drops out.

If `div` is omitted from argument list, a `div`>1 is calculated.

## Return codes

* `eu[1]==1` for existence
* `eu[2]==1` for uniqueness
* `eu[1]==-1` for existence only with not-s.c. z;
* `eu==[-2,-2]` for coincident zeros.

### Notes

We constrain Julia to use the complex version of the `schurfact` routine regardless of the
types of `Γ0` and `Γ1`, to match the behavior of Matlab.  Matlab always uses the complex version
of the Schur decomposition, even if the inputs are real numbers.
"""
function gensysdt(Γ0, Γ1, c, Ψ, Π, args...)
    F = schurfact!(complex(Γ0), complex(Γ1))
    gensysdt(F, c, Ψ, Π, args...)
end

function gensysdt(F::Base.LinAlg.GeneralizedSchur, c, Ψ, Π)
    gensysdt(F, c, Ψ, Π, new_div(F))
end

const ϵ = 1e-10 # sqrt(eps()) * 10
const _falvesCHECK = true

# Method that does the real work. Work directly on the decomposition F
function gensysdt(F::Base.LinAlg.GeneralizedSchur, c, Ψ, Π, div, )
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
        bigev, ueta, deta, veta = decomposition_svd(q2xΠ)
        _    , teta, seta, weta = decomposition_svd(q2xΨ)
    end

    # branch below to handle case of no stable roots
    if nunstab == n
        q1xΠ = zeros(0, neta)
        bigev1 = 0
        ueta1 = zeros(0, 0)
        veta1 = zeros(neta, 0)
        deta1 = zeros(0, 0)
    else
        bigev1, ueta1, deta1, veta1 = decomposition_svd(q1xΠ)
    end

    # Note that existence and uniqueness are not just matters of comparing
    # numbers of roots and numbers of endogenous errors.  These counts are
    # reported below because usually they point to the source of the problem.

    # REVIEW : picked from Reiter/Lubik
    if (_falvesCHECK)

        ## Existence (1)
        loose = teta - A_mul_Bc(ueta, ueta) * teta
        loosesvd = svdfact(loose)
        nloose = sum(abs(loosesvd[:S]) .> ϵ * n)
        iexist1  = (nloose == 0)

        ## Existence (2)
        rq2   = rank(q2xΠ, ϵ);
        rq2q2 = rank([q2xΨ q2xΠ], ϵ);
        iexist2  = rq2 == rq2q2;

        iexist = iexist1 && iexist2

        eu[1] = Int(iexist)

        # Uniqueness  -- eq 48, 49 on Sims
        if isempty(veta1)
            iunique = true
        else
            loose = veta1 - A_mul_Bc(veta, veta) * veta1
            loosesvd = svdfact(loose)
            nloose = sum(abs(loosesvd[:S]) .> ϵ * n)
            iunique = (nloose == 0)
        end

        eu[2] = Int(iunique)

        @printf("   Number of unstable eigenvalues : %d \n", nunstab)
        @printf("   q2 = (%d,%d); Ψ = (%d,%d); Π = (%d,%d)\n\n",
        size(qt2,2),size(qt2,1), size(Ψ,1),size(Ψ,2), size(Π,1),size(Π,2));

        if !all(eu.==1)
            iexist  || info("NO existence")
            iunique || info("Indeterminacy. $(nloose) loose endogeneous errors")
        end

        # eu[1] = length(bigev) <= neta                              # falves : from Miao textbook
        # maxabs( (eye(nunstab) - q2xΠ*((q2xΠ'*q2xΠ)\q2xΠ'))*q2xΨ )  # falves : residuals of regressing are zero
        # eu[1] = Int( rank(q2xΠ) == rank(hcat(q2xΨ,q2xΠ)) )         # falves : from checkeu Reiter

        ###  NOTE:  column rank = row rank  ###

    # original
    else
        # Existence
        iexist = length(bigev) >= nunstab

        if iexist
            eu[1] = 1
        else
            info("NO existence. $(nunstab) eigenvalues")
            @printf("Q = (%d,%d); q2 = (%d,%d); psi = (%d,%d); pi = (%d,%d)\n\n",
            size(qt,2),size(qt,1), size(qt2,2),size(qt2,1), size(Ψ,1),size(Ψ,2), size(Π,1),size(Π,2));
        end

        # Uniqueness
        if isempty(veta1)
            unique = true
        else
            loose = veta1 - A_mul_Bc(veta, veta) * veta1
            loosesvd = svdfact(loose)
            nloose = sum(abs(loosesvd[:S]) .> ϵ * n)
            unique = (nloose == 0)
        end

        if unique
            eu[2] = 1
        else
            info("Indeterminacy. $(nloose) loose endogeneous errors")
        end
    end
    Phi = ( ueta1 * A_mul_Bc(deta1, veta1) ) * ( veta * (deta \ ueta') )
    # ( ueta * (deta \ veta') * veta1 * A_mul_Bc(deta1, ueta1))'
    tmat = hcat( eye(n - nunstab), -Phi ) # Sims : I - Φ

    G0 = vcat(tmat * a, hcat(zeros(nunstab, n - nunstab), eye(nunstab)))
    G1 = vcat(tmat * b, zeros(nunstab, n))

    # G0 is always non-singular because by construction there are no zeros on
    # the diagonal of a(1:n-nunstab,1:n-nunstab), which forms G0's ul corner.
    G0I = inv(G0)
    G1 = G0I * G1
    usix = (n - nunstab + 1):n
    Busix = b[usix,usix]
    Ausix = a[usix,usix]
    C = G0I * vcat(tmat * Ac_mul_B(qt, c), (Ausix - Busix) \ Ac_mul_B(qt2, c))
    impact = G0I * vcat(tmat * Ac_mul_B(qt, Ψ), zeros(nunstab, size(Ψ, 2)))
    fmat = Busix \ Ausix
    fwt = -Busix \ Ac_mul_B(qt2, Ψ)
    ywt = G0I[:, usix]

    loose = G0I * vcat(q1xΠ * (eye(neta) - A_mul_Bc(veta, veta)), zeros(nunstab, neta))

    G1 = real(z * A_mul_Bc(G1, z))
    C = real(z * C)
    impact = real(z * impact)
    loose = real(z * loose)
    ywt = z * ywt

    return G1, C, impact, fmat, fwt, ywt, gev, eu, loose
end


function new_div(F::Base.LinAlg.GeneralizedSchur)
    a, b = F[:S], F[:T]
    n = size(a, 1)
    div = 1.01
    for i in 1:n
        if abs(a[i, i]) > 0
            divhat = abs(b[i, i] / a[i, i])
            if (1 + ϵ < divhat) && (divhat <= div)
                div = 0.5 * (1 + divhat)
            end
        end
    end
    return div
end

function decomposition_svd(A)
    Asvd = svdfact(A)               # REVIEW : don't change matrix A
    bigev = find(Asvd[:S] .> ϵ)
    Au = Asvd[:U][:, bigev]
    Ad = diagm(Asvd[:S][bigev])
    Av = Asvd[:V][:, bigev]
    return bigev, Au, Ad, Av

end
