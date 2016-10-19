
%%% Description
%       Holds all equilibrium conditions
%
%%% INPUTS
%       (1) X : Variables at t, t-1, expectational shocks ans regular shocks
%       (2) Env structure
%
%%% OUTPUT
%       (1) resid       : residual on all equilibrium conditions
%       (2) updated Env
%
function [resid, EnvOptional] = equHistogram(X, Env)

    global MP;

    xcurr = X(Env.iVarX.x);
    xlag  = X(Env.iVarX.xlag);
    eta   = X(Env.iVarX.eta);
    eps   = X(Env.iVarX.eps);

    % x2par(x,Env) OUTPUT
    %       vHistogram    : distribution of wealth
    %       Par  : policy rule parameters
    %       dz   : exog values
    %       KAggr: aggregate capital
    %       Env  : updated Env

    [vHistogram   , mSavingsPar   , dz   , KAggr, EnvOut]  = unpackHistogram_X(xcurr, Env);
    [vHistogramlag, mSavingsParlag, dzlag, KAggrlag]       = unpackHistogram_X(xlag , Env);

    %== Recover PRICES ==%
    R      = 1 + netintr(KAggr, 1+dz);
    wage   = wagefunc(KAggr, 1+dz);

    Rlag      = 1 + netintr(KAggrlag, 1 + dzlag);
    wagelag   = wagefunc(KAggrlag, 1 + dzlag);

    % ******************************************************************
    %   EQUATIONS
    % ==================================================================

    % Exogenoous Shocks
    resZ = [dz - MP.rhoTFP * dzlag - MP.sigTFP * eps(1)];       % REVIEW : does MP.sigTFP enters here?

    %-----------------------------------------------%
    %%% Dynamics of distribution of wealth        %%%
    %-----------------------------------------------%
    % tic
    %%%  NOTE:  using mSavingsParlag on this version of code   %%%
    %%%         vHistogram is the distribution over k_i and θ_j today

    %== Transition matrix ==%
    Pi = forwardmat(0, mSavingsParlag);

    %== Iterate one step ==%
    vHistogramLoM = forward(vHistogramlag, Pi);

    % toc
    %== Compute Transition matrix (iFullMat=1) ==%
    % Pi_ifull = forwardmat(1, mSavingsParlag);
    % D2_ifull = forward(vHistogramlag,Pi_ifull);
    % toc
    %== CHECK ==%
    % all( all(full(vHistogramLoM == D2_ifull)) )

    % LAW of motion for moments
    resD = distr2xdistr(vHistogramLoM, Env.vHistogramSS) - distr2xdistr(vHistogram, Env.vHistogramSS);

    % household policy rules
    resC = eulerres(mSavingsParlag, mSavingsPar, Rlag, R, wagelag, wage) + eta;

    if (MP.iIncludeKAggr)

        resK = KAggr - expect_k(vHistogram);
    else

        resK = [];
    end

    resid = [resD;          % Distribution
             resZ;          % Exogenous shocks
             resK;          % Aggregate capital
             resC];         % Household policies


    %== UPDATE Environment ==%
    if nargout >1

        EnvOptional = EnvOut;
    end
end
