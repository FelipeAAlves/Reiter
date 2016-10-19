% ******************************************************************
%   LINEARIZED Model (Reiter)
% ==================================================================

%%%  WARN:  CHECK which variables are eval at zero and which are   %%%
%%%         at their stst level. CONFIRM this inside equ.m         %%%
if (MP.iIncludeKAggr)
    X = [zeros(MP.nHistogramTotal-1,1); zeros(MP.nz,1); Kequ; Par(:)];
else
    X = [zeros(MP.nHistogramTotal-1,1);zeros(MP.nz,1);Par(:)];
end
MP.ntotal = length(X);

%== Prepate X ==%
[X2, varindx] = setix(X, MP.nSavingsPar*MP.neps, MP.nz);

%== UPDATE Env with Kequ and varindx ==%
Env.iVarX = varindx;

%== CHECK stst ==%
fprintf('STEP [1]: Checking stst K = %.6f\n', Kequ)
[residStst,Env] = equ(X2,Env);                  % UPDATE Env variable
if (max(abs(residStst))>1e-6)
    error('wrong stst');
end

%== Linearizing the MODEL ==%
disp('')
disp('STEP [2]: Linearizing model equations')
njac = 20*ones(size(Env.iVarX.x));
njac(Env.iVarAggr)=500;

Env = dojacob(@equ,X2,Env,njac,Env);

%== Solving the MODEL ==%
disp('')
disp('STEP [3]: Solving the model')
G1 = []; impact = [];abseig = [];
[G1,impact,eu,abseig] = solveexact(Env);
fprintf('[eu1, eu2] = %d, %d \n', eu(1),eu(2))
if (~all(eu==1))  %signal failure:
    G1 = [];
    warning('sims failed');
end

%====================================%
%%%      IRF and Simulation        %%%
%====================================%

%== Create observation matrix for all aggregate variables ==%
iA  = Env.iVarAggr;
niA = length(iA);

Hagg = sparse(1:niA, iA, ones(1,niA), niA, MP.ntotal);

irf_agg = ir_sims(G1, impact, MP.IRF_Length-1, Hagg,MP.Sigma);

% disp('Computing variances through simulation')
[ser, shocks] = simulateSystem(G1, impact, Hagg, MP.sim_T, MP.Sigma, MP.seed);

%== Recover Aggregate series ==%
z   = ser(:,1);
tau = ser(:,2);
K = ser(:,3) + Env.Kstst;


% ******************************************************************
%   FIGURES
% ==================================================================

figure
plot(ser(1:100,:), 'LineWidth', 3)
legend('Tech.', 'Taxes', 'Capital', 'Location','Best')

% .....................................................................................
