%%% Description:
%       Unpack the vector X (variable values) into its
%       components in Polynomial case
%
%%% INPUTS:
%       X vector of values
%
%%% OUTPUT:
%       (1) mMoments     : moments of wealth distribution
%       (2) mDensityCoeff  : Parameters of the density function
%       (3) dz
%       (4) KAggr
%       (5) mSavingsPar
%
%
%
function [mMoments, mDensity, dz, KAggr, mSavingsPar, Env ] = unpackPolynomial_X(X , Env)

    global MP;

    %== Sizes ==%
    nX      = length(X);
    neps    = MP.neps;

    nMom    = MP.nMoments;
    nMomTt  = MP.nMomentsTotal;

    nDensExp   = MP.nDensityCoeff-1;   % take out 1 bc X does not include g₀
    nDensExpTt = nDensExp * neps;

    nAggr   = MP.nAggr;
    nc      = MP.nSavingsPar;

    %== Recover Moments ==%
    iVarMom  = 1:nMomTt;
    mMoments = reshape( X(iVarMom), nMom, neps );

    %== Recover Density Parameters ==%
    tmp_n0   = nMomTt + nDensExpTt;
    iVarDens = nMomTt+1:tmp_n0;

    [mDensity, ~] = XtoDensity( mMoments, X(iVarDens) );

    %== Recover Aggregates ==%
    tmp_n1 = nMomTt + nDensExpTt + nAggr;

    iVarAggr  = [tmp_n0+1:tmp_n1];
    iVarAggr_exog = tmp_n0 + find(MP.AggrInd.exog);
    iVarAggr_endo = tmp_n0 + find(MP.AggrInd.endo);

    iVarState     = [1:tmp_n0 tmp_n0+find(MP.AggrInd.state)];


    % Shock
    dz        = X( iVarAggr_exog );

    % Aggregate Capital
    if (MP.iIncludeKAggr)
        KAggr = X(iVarAggr_endo);
    end

    %== Savings Policy ==%
    tmp_n0 = tmp_n1;
    iVarHousehold = (tmp_n0 + 1):nX;

    mSavingsPar           = reshape( X(iVarHousehold), nc, neps );

    if nargin==2 && nargout>5
        % Moments
        Env.iVarMom         = iVarMom;

        % Density params
        Env.iVarDens        = iVarDens;

        % Aggr
        Env.iVarAggr        = iVarAggr;
        Env.iVarAggr_exog   = iVarAggr_exog;
        Env.iVarAggr_endo   = iVarAggr_endo;
        Env.iVarState       = iVarState;

        % Policies
        Env.iVarHousehold   = iVarHousehold;
    end
