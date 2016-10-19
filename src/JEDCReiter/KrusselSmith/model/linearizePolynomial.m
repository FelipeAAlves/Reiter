% ******************************************************************
%   LINEARIZED Model (Winberry)
% ==================================================================

%%%  WARN:  CHECK which variables are eval at zero and which are   %%%
%%%         at their stst level. CONFIRM this inside equ.m         %%%

%== Extract only the Exp terms... Adjust g₀ inside equ to make sure ∫ to 1 ==%
mDensityExpCoeff = SS_Polynomial.mDensityCoeff(2:end,:);

if (MP.iIncludeKAggr)
    X = [ SS_Polynomial.mMoments(:); mDensityExpCoeff(:); zeros(MP.nz, 1); SS_Polynomial.K; SS_Polynomial.mSavingsPar(:) ];
else
    X = [ SS_Polynomial.mMoments(:); mDensityExpCoeff(:); zeros(MP.nz, 1); SS_Polynomial.mSavingsPar(:)];
end

Env.ntotal = length(X);

%== Prepate X ==%
% setiX(X, nEta, nEps)
[X2, Env.iVarX] = setiX(X, MP.nSavingsPar * MP.neps, MP.nz);        % 1st instance of Env


fprintf('STEP [1]: Checking stst K = %.6f\n', SS_Polynomial.K)

%== CHECK stst ==%
[residStst, Env] = equPolynomial(X2,Env);                  % UPDATE Env variable
if (max(abs(residStst))>1e-6)
    error('wrong stst');
end

disp('')
disp('STEP [2]: Linearizing model equations')

%== Linearizing the MODEL ==%
njac = 20*ones(size(Env.iVarX.x));

njac(Env.iVarMom)  = 1000/(MP.nMomentsTotal);
njac(Env.iVarDens) = 1000/(MP.nDensityCoeffTotal);
njac(Env.iVarAggr) = 500;

Env = dojacob(@equPolynomial,X2,Env,njac,Env);

disp('')
disp('STEP [3]: Solving the model')

%== Solving the MODEL ==%
G1 = []; impact = [];abseig = [];
[G1,impact,eu,abseig] = solveexact(Env);
fprintf('[eu1, eu2] = %d, %d \n', eu(1),eu(2))
if (~all(eu))  %signal failure:
    G1 = [];
    warning('sims failed');
end

%====================================%
%%%      IRF and Simulation        %%%
%====================================%

%== Construct H observation matrix ==%
nMom  = length(Env.iVarMom);
nDens = length(Env.iVarDens);
nAggr = length(Env.iVarAggr);

nTotalObs = nMom + nDens + nAggr;
iHobs     = [Env.iVarMom Env.iVarDens Env.iVarAggr];

Hobs = sparse(1:nTotalObs, iHobs, ones(1,nTotalObs), nTotalObs, Env.ntotal);

%== Series Indicators ==%
iSeriesPoly.Mom  = logical([ones(1,nMom) zeros(1,nTotalObs-nMom)]);
tmp_n = nMom;
iSeriesPoly.Dens = logical([zeros(1,tmp_n) ones(1,nDens) zeros(1,nTotalObs-tmp_n-nDens)]);
tmp_n = tmp_n + nDens;
iSeriesPoly.Aggr = logical([zeros(1,tmp_n) ones(1,nAggr)]);

% Impulse Response
% [series, shocks] = SimulateSystem(G1, impact, MP.IRF_Length, MP.Sigma, 'irf');

% Simulation
[series, shocks] = SimulateSystem(G1, impact, MP.sim_T, MP.Sigma, 'simulation', MP.seed);

seriesPolynomial = (Hobs * series)';
clearvars series
