% INITIALIZE Parameter Values

%%% Description:
%       Creates a global struc MP to hold all relevant Model NOTErmation
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

MP.ydist         = [1 - aggEmployment;aggEmployment];                            % distribution
MP.aggEmployment = dot(MP.ydist, MP.ypp);                                        % aggregate Employment

MP.mu = .15;                                                                     % unemployment benefits
MP.tau = MP.mu * (1 - aggEmployment) / aggEmployment;                            % flat tax

% wage already contain taxes
% MP.iWageTax = 1; DEPRECATED

% ******************************************************************
%   MODEL PARAMETERS:
% ==================================================================
% borrowing constraint:
MP.aabar = 0;

% Economy Parametrization
MP.alpha = 0.36;
MP.beta  = 0.96;
MP.delta = 0.1;
MP.gam   = 1;

% Representative AGENT Stst
%%% OPTION 01 : Normalize stst capital to 1 %%%
%%%  WARN:  Having problem afterwards when I need to compute map m --> g %%%

% R        = 1/MP.beta;
% MP.A = ( (R-1+MP.delta) / MP.alpha ) * MP.aggEmployment ^ (MP.alpha - 1);
% KLratio  = inv( MP.aggEmployment );
% assert(abs(R-1-netintr(1.0,1))<1e-8);

%%% OPTION 02 : Traditional Problem %%%
MP.A = 1;
MP.KRepSS = ( ( MP.alpha * (MP.aggEmployment ^ (1 - MP.alpha)) ) / ( (1 / MP.beta) - (1 - MP.delta) ) ) ^ (1 / (1 - MP.alpha));

% exog shocks
MP.nz = 1;

% correlation shocks:
MP.sigTFP  = 0.014;
MP.rhoTFP  = 0.859;

MP.Sigma = [1];                 % IMPORTANT: Set it to 1... otherwise I am counting twice...
MP.rhozz = diag([MP.rhoTFP]);

% ******************************************************************
%   IMPLEMENTATION CHOICES
% ==================================================================

% LIMITS on Capital capital:
MP.AssetsMin = MP.aabar;
MP.AssetsMax = 3 * MP.KRepSS;                               % original = 3, CHANGE line 127 if increased to have 1.5 MP.KRepSS
                                                            %               CHANGE nAssetsQuad to compute integral
%MP.AssetsMin = MP.aabar;  MP.AssetsMax = 5;

% number of grid points for each shock:
MP.nSavingsPar    = 120;                    % policy function; default = 80-160
MP.nHistogram     = 75;                     % finer grid     ; default = 75-100
                                            %       bigger values perform even worst during linearization
                                            %       bc what changes is only the weigth bet adjancent points but not the POINT where I go

% APPROXIMATION of density
MP.nMoments     = 4;                                                        % original: 3
MP.nDensityCoeff = MP.nMoments+1;                                           % density coefficients
MP.nAssetsQuad  = 12;                                                       % original: 8,
                                                                            %           12 works for MP.AssetsMax = 3 * MP.KRepSS
                                                                            %           16 works for MP.AssetsMax = 4 * MP.KRepSS
                                                                            %           20 works for MP.AssetsMax = 5 * MP.KRepSS
MP.nMomentsTotal = MP.neps * MP.nMoments;                                           % total moments
MP.nDensityCoeffTotal = MP.neps * MP.nDensityCoeff;                                 % total coeff


% use KAggr as extra variable
MP.iIncludeKAggr = 1;

MP.nAggr = MP.iIncludeKAggr + MP.nz;

%== CREATE indices ==%
% MP.aggInd.dec    = logical([zeros(1,MP.nagg-n) ones(1,n)]);
MP.AggrInd.state  = logical([1 1]);
MP.AggrInd.exog   = logical([1 0]);
MP.AggrInd.endo   = logical([0 1]);
% MP.aggInd.static = ~logical(MP.aggInd.state + MP.aggInd.dec);

% ******************************************************************
%   FINE grid for wealth histogram:
% ==================================================================
MP.nHistogramTotal = MP.nHistogram*MP.neps;  %makes only sense for neps==1; otherwise_:DSF!

% Accumulate points in the beginning
% [MP.AssetsGridFine, MP.logshift] = makeknotd(MP.AssetsMin,MP.AssetsMax,MP.nHistogram);
% [temp, MP.logshift] = makeknotd(MP.AssetsMin,MP.AssetsMax,MP.nHistogram);
MP.logshift = 1.0;

% Wiberry GRID
MP.AssetsGridFine = linspace(MP.AssetsMin,MP.AssetsMax,MP.nHistogram)';

% ******************************************************************
%   Quadrature grid, for computing integrals
% ==================================================================

% compute quad points and weight
[vAssetsGridQuadZeros, MP.QuadWeights] = computeGaussLegendreQuadrature(MP.nAssetsQuad);

% scale to the assetGrid
MP.AssetsGridQuad = scaleUp(vAssetsGridQuadZeros, MP.AssetsMin + 1e-1, MP.AssetsMax);


% ******************************************************************
%   Knot points for EULER Collocation:
% ==================================================================

MP.xmin = 0.000;
MP.xmax = MP.AssetsMax * 3/8; % beyond that point, linear extrapolation is fine!

n1        = ceil(MP.nSavingsPar/3);
x1        = linspace(MP.xmin,0.25,n1+1)';
x2        = logspaceshift(0.25, MP.xmax, MP.nSavingsPar-n1, MP.logshift)';
EulerGrid = [x1(1:end-1);x2];
MP.knotXi = EulerGrid(2:end);   %take out the zero

%  MP.knotXi = logspaceshift(MP.xmin,MP.xmax,MP.nSavingsPar,MP.logshift)';

% initial GUESS parameter
Par0 = [0.0; 0.8*MP.knotXi];
MP.SavingsParstart = repmat(Par0, MP.neps,1);

% use polynomials in logs, for moments:
MP.momtype = 'pl';

% ******************************************************************
%   Options for the solution algorthim
% ==================================================================
MP.IRF_Length = 40;
MP.sim_T      = 1000;
MP.seed       = 4239074;
