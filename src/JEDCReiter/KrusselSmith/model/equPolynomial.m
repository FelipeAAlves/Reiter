
%%% Description
%       Holds all equilibrium conditions
%
%%% INPUTS
%       (1) X : variables t, t-1, exp shocks ans reg shocks
%       (2) Env structure
%
%%% OUTPUT
%       (1) resid
%       (2) updated Env
function [resid, EnvOptional] = equPolynomial(X, Env)

    global MP;

    xcurr = X(Env.iVarX.x);
    xlag  = X(Env.iVarX.xlag);
    eta   = X(Env.iVarX.eta);
    eps   = X(Env.iVarX.eps);

    %== Unpack from X ==%
    [mMoments   ,mDensity, dz, KAggr, mSavingsPar, EnvOut ]    = unpackPolynomial_X(xcurr, Env);
    [mMomentslag,mDensitylag, dzlag, KAggrlag, mSavingsParlag] = unpackPolynomial_X(xlag);

    % recover Prices
    R      = 1 + netintr(KAggr, 1+dz);
    wage   = wagefunc(KAggr, 1+dz);
    Rlag      = 1 + netintr(KAggrlag, 1 + dzlag);
    wagelag   = wagefunc(KAggrlag, 1 + dzlag);

    % ******************************************************************
    %   EQUATIONS
    % ==================================================================

    % Exogenoous Shocks
    resZ = [dz - MP.rhoTFP * dzlag - eps(1)];

    %------------------------------------------------------------------------%
    %%% Dynamics of Moments and Consistency of Moments × Parameters        %%%
    %------------------------------------------------------------------------%

    % Sizes
    nAssetsQuad = MP.nAssetsQuad;
    nMoments = MP.nMoments;
    neps     = MP.neps;

    % Grids
    AssetsGridQuad = MP.AssetsGridQuad;

    %== Construct PREVIOUS period policy eval Quad nodes ==%
    mAssetsPrime = initsize(zeros(nAssetsQuad, neps), mSavingsParlag);
    for ieps = 1 : neps

        %== Saving policy ==%
        S = savingspline( mSavingsParlag(:,ieps) );

        %== Evaluate at the Quad nodes ==%
        mAssetsPrime(:,ieps) = interp_savspline(S, AssetsGridQuad);

    end

    %===================================%
    %%% Consistency CONDITIONS        %%%
    %===================================%
    % aGridMoments  = initsize( zeros(nAssetsQuad, nMoments), mMomentslag);
    % mDensitylag2  = initsize( zeros(nAssetsQuad, neps), mDensityCoefflag, mMomentslag);
    mMomentslag2  = initsize( zeros(size(mMomentslag)), mDensitylag );

    for ieps =1:neps

        %======================================================%
        %%% Density recovered together with variables        %%%
        %======================================================%

        % %== Set Up aGrid ==%
        % aGridMoments(:,1) = ( AssetsGridQuad - mMomentslag(1,ieps) );
        %
        % for iMoment = 2 : nMoments
        %     aGridMoments(:,iMoment) = ...
        %                 ( AssetsGridQuad - mMomentslag(1,ieps) ) .^ iMoment - mMomentslag(iMoment,ieps);
        % end
        %
        % %== Recover Density ==%
        % mDensitylag2(:,ieps) = mDensityCoefflag(1,ieps) * exp( aGridMoments * mDensityCoefflag(2:end,ieps) );

        %== Moments ==%
        mMomentslag2(1,ieps) = MP.QuadWeights' * ( AssetsGridQuad .* mDensitylag(:,ieps) );

        for iMoment = 2 : nMoments
            mMomentslag2(iMoment,ieps) = ...
                    MP.QuadWeights' * ( ...
                    ( AssetsGridQuad - mMomentslag(1,ieps) ) .^ iMoment .* ...
                    mDensitylag(:,ieps) );
        end
    end

    resDens = mMomentslag(:) - mMomentslag2(:);

    %=====================================%
    %%% LAW OF MOTION of Moments        %%%
    %=====================================%

    mMoments2 = initsize( zeros(size(mMoments)), mAssetsPrime, mDensitylag, mMoments );
%     mMoments2 = initsize( zeros(size(mMoments)), mAssetsPrime, mDensitylag);

    for ieps = 1 : neps % WARN THIS period shock

        %%%  NOTE:  computing next period first moments     %%%
        %%% m' = ∑ π( ̃ϵ | ϵ ) ∫ a'( ̃ϵ , a; z, m ) g_ϵ(a) da

        LawofMotionAux = 0;
        for jeps = 1 : neps % WARN PREVIOUS period shock

            LawofMotionAux = LawofMotionAux + ...
                    MP.ydist(jeps) * MP.transpp(jeps, ieps) * ...
                    MP.QuadWeights' * ( mAssetsPrime(:,jeps) .* mDensitylag(:, jeps) );
        end
        mMoments2(1,ieps) = LawofMotionAux / MP.ydist(ieps);


        % Compute higher order moments (centered)
        for iMoment = 2 : nMoments

            %%%  NOTE:  computing next period ith moment                    %%%
            %%%         Moment inside is with respect to ieps               %%%
            %%%
            %%%     mᵢ'(ϵ) = ∑ π( ̃ϵ | ϵ ) ∫ [ a'( ̃ϵ , a; z, m ) - m��?(ϵ) ]^i g_̃ϵ (a) da
            %%%

            LawofMotionAux = 0;
            for jeps = 1 : neps % WARN THIS period shock

                LawofMotionAux = LawofMotionAux + ...
                    MP.ydist(jeps) * MP.transpp(jeps, ieps) * ...
                    MP.QuadWeights' * ( ...
                    ( ( mAssetsPrime(:, jeps) - mMoments(1,ieps) ) .^ iMoment ) .* ...
                    mDensitylag(:,jeps) );
            end

            mMoments2(iMoment, ieps) = LawofMotionAux / MP.ydist(ieps);

        end

    end

    resMom = mMoments(:) - mMoments2(:);

    % household policy rules
    resConsump = eulerres(mSavingsParlag, mSavingsPar, Rlag, R, wagelag, wage) + eta;

    if (MP.iIncludeKAggr)

        resK = KAggr - mMoments(1,:) * MP.ydist;
    else

        resK = [];
    end

    resid = [resDens;               % Consistency of density × Moments
             resMom;                % Law of Motion of Moments
             resZ;                  % Exogenous shocks
             resK;                  % Aggregate capital
             resConsump];           % Household policies


    %== UPDATE Environment ==%
    if nargout >1

        EnvOptional = EnvOut;
    end
end
