
% INITIALIZE Parameter Values

%%% Description:
%       Creates a struc MP to hold all relevant NOTErmation
%
global MP;

% ******************************************************************
%   DISCRETIZATION OF INCOME PROCESS
% ==================================================================
MP.neps = 2;                             % permanent productivity states

MP.ypp = [0 1.0];
uDuration = 1;                              %number of period unemployed before employed
aggEmployment  = .93;

% Transition prob
MP.transpp = [uDuration / (1 + uDuration), 1 - (uDuration / (1 + uDuration));...
((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration))),...
1 - ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration)))...
];

MP.ydist         = [1 - aggEmployment;aggEmployment];                   % distribution
MP.aggEmployment = dot(MP.ydist, MP.ypp);                       % aggregate Employment

MP.mu = .15;                                                                  % unemployment benefits
MP.tau = MP.mu * (1 - aggEmployment) / aggEmployment;                            % flat tax

% ******************************************************************
%   MODEL PARAMETERS:
% ==================================================================

% borrowing constraint:
MP.AssetsMin = 0;

% wage already contain taxes
MP.iWageTax = 1;

% Model parameters:
MP.alpha = 0.36;
MP.beta  = 0.96;
MP.delta = 0.1;
MP.gam   = 1;

% Representative AGENT Stst
% take into account that Y does not add up exactly to 1 because of truncation
R        = 1/MP.beta;
MP.A = ( (R-1+MP.delta) / MP.alpha ) * MP.aggEmployment ^ (MP.alpha - 1);
KLratio  = inv( MP.aggEmployment );
assert(abs(R-1-netintr(1.0,1))<1e-8);

% Upper limit capital:
MP.AssetsMax = 10;

% exog shocks
MP.nz = 1;

% correlation shocks:
MP.sigTFP  = 0.014;
MP.rhoTFP  = 0.859;

MP.Sigma = 0;
MP.rhozz = diag([MP.rhoTFP]);

% use KAggr as extra variable
MP.iIncludeKAggr = 1;

MP.nagg = MP.iIncludeKAggr + MP.nz;

%== CREATE indices ==%
% MP.aggInd.dec    = logical([zeros(1,MP.nagg-n) ones(1,n)]);
MP.aggInd.state  = logical([1 1 1]);
MP.aggInd.exog   = logical([1 1 0]);
MP.aggInd.endo   = logical([0 0 1]);
% MP.aggInd.static = ~logical(MP.aggInd.state + MP.aggInd.dec);

% number of grid points for each shock:
MP.nSavingsPar     = 100;            % policy function
MP.nHistogram = 500;            % finer grid

% APPROXIMATION of density
MP.nMoments     = 3;
MP.nAssetsQuad  = 8;
MP.nMomentsCoefficients = MP.neps * MP.nMoments;


% ******************************************************************
%   FINE grid for wealth histogram:
% ==================================================================
MP.nHistogramTotal = MP.nHistogram*MP.neps;  %makes only sense for neps==1; otherwise_:DSF!

% number knot points:
ndk = MP.nHistogram;
[MP.AssetsGridFine, MP.logshift] = makeknotd(MP.AssetsMin,MP.AssetsMax,ndk);


% ******************************************************************
%   Quadrature grid, for computing integrals
% ==================================================================

[AssetsGridQuadZeros, MP.QuadWeights] = computeGaussLegendreQuadrature(MP.nAssetsQuad);

% scale to the assetGrid
MP.AssetsGridQuad = scaleUp(vAssetsGridQuadratureZeros, MP.AssetsMin +1e-1,MP.AssetsMax);


% ******************************************************************
%   Knot points for EULER Collocation:
% ==================================================================

MP.xmin = 0.000;
MP.xmax = MP.AssetsMax*0.5; % beyond that point, linear extrapolation is fine!

n1 = ceil(MP.nSavingsPar/3);
x1 = linspace(MP.xmin,0.2,n1+1)';
x2 = logspaceshift(0.2,MP.xmax,MP.nSavingsPar-n1,MP.logshift)';
EulerGrid = [x1(1:end-1);x2];
MP.knotXi = EulerGrid(2:end);   %take out the zero

%  MP.knotXi = logspaceshift(MP.xmin,MP.xmax,MP.nSavingsPar,MP.logshift)';

% initial GUESS parameter
Par0 = [0.0; 0.95*MP.knotXi];
MP.Parstart = repmat(Par0, MP.neps,1);


% use polynomials in logs, for moments:
MP.momtype = 'pl';

% ******************************************************************
%   INDEXING
% ==================================================================


% ******************************************************************
%   Options for the solution algorthim
% ==================================================================
MP.IRF_Length = 40;
MP.sim_T = 1000;
MP.seed = 4239074;
