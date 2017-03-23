
function klein(∇f_x′, ∇f_y′, ∇f_x, ∇f_y, stake=1.0)

    Γ0 =  [ ∇f_x′ ∇f_y′ ];
    Γ1 = -[ ∇f_x ∇f_y ];
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

function findnzrows(A,eps=0.0)
  return find(sum(abs(A),2).>eps)
end

function findnzcols(A,eps=0.0)
  return find(sum(abs(A),1).>eps)
end

function nzcols(A,eps=0.0)
  return sum(sum(abs(A),1).>eps)
end
