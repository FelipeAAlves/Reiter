function [G1,impact,eu,abseig] = solveexact(Env)

    n   = size(Env.jac0,1);

    g0  = -Env.jac0;                  % current variables
    g1  =  Env.jac1;                  % lagged variables
    c   = zeros(size(g0,1),1);        % constants

    Psi = Env.jace;                   % Psi multiplies shocks in Sims paper
    Pi  = Env.jact;                   % expectational errors

    div = 1+1e-10;

    cput = cputime;
    [eu,eigvals,G1,impact] = checkeu(full(g0),full(g1),full(Psi),full(Pi),div);
    abseig = flipud( sort(abs(eigvals)) );

    if (any(eu==0))
        warning(sprintf('error in solveexact; eu = %d %d',eu(1),eu(2)));
    end
