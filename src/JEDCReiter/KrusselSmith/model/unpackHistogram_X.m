%%% Description:
%       Unpack the vector X (variable values) into its
%       components (distribution of wealth, policy rule parameters, aggregate
%       variables)
%
%%% INPUT:
%       X vector of values
%
%%% OUTPUT:
%       vHistogram    : distribution of wealth
%       mSavingsPar  : policy rule parameters
%       dz   : exog values
%       xaggr: aggregate variable values
%
%
function [vHistogram,mSavingsPar,dz,KAggr,Env] = unpackHistogram_X(X,Env)

    global MP;

    %== Sizes ==%
    nX   = length(X);
    nHistogramX   = MP.nHistogramTotal - 1;       % nStatesDistr = ndstst*neps

    nAggr = MP.nAggr;

    neps = MP.neps;
    nSavingsPar   = MP.nSavingsPar;

    %== Recover vHistogram ==%
    iVarHistogram = 1:nHistogramX;
    vHistogram = x2distr( X(iVarHistogram), Env.vHistogramSS );

    tmp_n = nHistogramX + nAggr;

    %== Set indexes in Env ==%
    iVarAggr  = [nHistogramX+1:tmp_n];

    iVarAggr_exog = nHistogramX + find(MP.AggrInd.exog);
    iVarAggr_endo = nHistogramX + find(MP.AggrInd.endo);
    iVarState     = [1:nHistogramX nHistogramX+find(MP.AggrInd.state)];

    iVarHousehold       = (tmp_n + 1):nX;

    % Shock
    dz        = X( iVarAggr_exog(1) );

    KAggr = X(iVarAggr_endo);

    %== recover mSavingsPar ==%
    mSavingsPar           = reshape( X(iVarHousehold), nSavingsPar, neps);

    if nargout>4

        %== Histogram ==%
        Env.iVarHistogram   = iVarHistogram;

        %== Aggregates ==%
        Env.iVarAggr        = iVarAggr;
        Env.iVarAggr_exog   = iVarAggr_exog;
        Env.iVarAggr_endo   = iVarAggr_endo;
        Env.iVarState       = iVarState;

        %== policies ==%
        Env.iVarHousehold   = iVarHousehold;
    end

    % if (MP.iIncludeKAggr)
    %
    %     % Env.iVarStatic  = nd + find(MP.AggrInd.static);
    %     KAggr = X(Env.iVarAggr_endo);
    %
    % else
    %
    %     % Env.iVarStatic = [];
    %     % if (MP.iKS)
    %     %     % KAggr = vHistogram(1)*Env.norm_h1;
    %     %     KAggr = Env.mom2K(1:nd)'*X(1:nd);
    %     % else
    %     %     KAggr = make_h(1)' * X(1:nd) + Env.Kstst;
    %     % end
    %
    % end
