% ******************************************************************
%   LINEARIZED Model (Reiter)
% ==================================================================

%%%  WARN:  CHECK which variables are eval at zero and which are   %%%
%%%         at their stst level. CONFIRM this inside equHistogram.m         %%%
if (MP.iIncludeKAggr)
    X = [zeros(MP.nHistogramTotal-1,1); zeros(MP.nz,1); SS_Histogram.K; SS_Histogram.mSavingsPar(:)];
else
    X = [zeros(MP.nHistogramTotal-1,1); zeros(MP.nz,1); SS_Histogram.mSavingsPar(:)];
end

Env.ntotal = length(X);

%== Prepate X and Env ==%
[X2, Env.iVarX] = setiX(X, MP.nSavingsPar * MP.neps, MP.nz);  %1st instace of Env Structure

%== ADD vHistogram to Env variable - WILL be used inside the equHistogram ==%
Env.vHistogramSS = SS_Histogram.vHistogram;

fprintf('STEP [1]: Checking stst K = %.6f\n', SS_Histogram.K)

%== CHECK StSt ==%
[residStst, Env] = equHistogram(X2, Env);                  % UPDATE Env variable

if (max(abs(residStst))>1e-6)
    error('wrong stst');
end

disp('')
disp('STEP [2]: Linearizing model equations')

%== Linearizing the MODEL ==%
njac = 20 * ones( size(Env.iVarX.x) );
njac( Env.iVarAggr ) = 500;

Env = dojacob( @equHistogram, X2, Env, njac, Env );

disp('')
disp('STEP [3]: Solving the model')

%== Solving the MODEL ==%
G1 = []; impact = [];abseig = [];
[G1,impact,eu,abseig] = solveexact(Env);
fprintf('[eu1, eu2] = %d, %d \n', eu(1),eu(2))
if (~all(eu==1))  %signal failure:
    G1 = [];
    warning('sims failed');
end

% ******************************************************************
%   IRF and SIMULATION
% ==================================================================

%== Create observation matrix for all aggregate variables ==%
nHistogramX = length(Env.iVarHistogram);
nAggr       = length(Env.iVarAggr);
nSavingsPar = length(Env.iVarHousehold);

nTotalObs = nHistogramX + nAggr + nSavingsPar;
iHobs     = [Env.iVarHistogram Env.iVarAggr Env.iVarHousehold];

Hobs = sparse(1:nTotalObs, iHobs, ones(1,nTotalObs), nTotalObs, Env.ntotal);

%== Series Indicators ==%
iSeriesHist.HistogramX = logical([ones(1,nHistogramX) zeros(1,nTotalObs-nHistogramX)]);

tmp_n = nHistogramX;
iSeriesHist.Aggr = logical([zeros(1,tmp_n) ones(1,nAggr) zeros(1,nTotalObs - tmp_n - nAggr)]);

tmp_n = tmp_n + nAggr;
iSeriesHist.SavingsPar = logical([zeros(1,tmp_n) ones(1,nTotalObs - tmp_n)]);

% HAggr = sparse(1:nAggr, Aggr, ones(1,nAggr), nAggr, Env.ntotal);

%---------------------%
%%% IRF        %%%
%---------------------%
[series, shocks] = SimulateSystem(G1, impact, MP.IRF_Length, 5*MP.Sigma, 'irf');

irfHist = (Hobs * series{1})';
clearvars series

%---------------------%
%%% SIMULATION        %%%
%---------------------%
[series, shocks] = SimulateSystem(G1, impact, MP.sim_T, MP.Sigma, 'simulation', MP.seed);

simulHist = (Hobs * series)';
clearvars series Env
