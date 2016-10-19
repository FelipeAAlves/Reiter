%%% Description:
%       Computes residual of market-clearing condition, the density as an approximation
%       for the distribution;
%
%%% INPUT
%       (1) capital: candidate aggregate capital stock
%       (2) mMoments:  (nMomentsn × neps ) intial guess of moments of distribution
%       (3) aGridMoments:  (nAssetsQuad x nMoments × neps) grid of centralized moments that enter on pdf, corresponding to the nMoments
%
%       (4) (optional) aggStatus: 1 only when solving for other quantities (not
%           just residual) in the full dynamic version
%
%%% OUTPUT
%   (1) residual: residual of market clearing condition

function [resid, mSavingsParOptional, mDensityCoeffOptional, mMomentsOptional, SS ] = ststPolynomialResid(K, mMoments, aGridMoments, print)


    if nargin == 3
        print = 0;
    end

    %%%  OPTION  %%%
    toldensity  = 1e-6;                     % WARN : anything smaller will activate the assert
    tolmoments  = 1e-7;                     % WARN: making it smaller than 1e-7 is not recommended
                                            %       since fminunc itself doesn't have that much precision
    maxiter     = size(mMoments,1) * 500;

    global MP;

    %== Sizes ==%
    nc   = MP.nSavingsPar;                       % grid on savings policy
    neps = MP.neps;                     % idio show
    nMoments    = MP.nMoments;
    nAssetsQuad = MP.nAssetsQuad;

    if print
        fprintf(1,' OUTER loop: \n')
        fprintf(1,'  K = %.3f \n', K);
    end
    %%%  NOTE:  given a GUESS for Kequ, recover prices           %%%
    %%%         and solve for consumption pol by collocation     %%%

    % recover prices
    R      = 1 + netintr(K,1);
    wage   = wagefunc(K, 1);

    %=======================================================%
    %%% Compute Policy function for set of prices         %%%
    %=======================================================%
    if print
        disp(['   INNER LOOP (policy)'])
    end

    %== Always restart at same GUESS ==%
    [vSavingsPar, check] = broydn(@eulerres_stst,MP.SavingsParstart,[1e-11,1,print],R,wage);
    if (check~=0)
        warning('broydn not converged');
    end

    mSavingsPar = reshape(vSavingsPar, nc, neps);

    %=================================================================%
    %%% Compute Stationary Distribution from decision rules         %%%
    %=================================================================%
    if print
        disp(['   FINDING Stationary Distribution'])
        % fprintf('\n   ')
    end

    %== Compute savings Policy on Quad grid ==%
    mAssetsPrime = zeros(nAssetsQuad, neps);
    for ieps = 1 : neps

        %== Saving policy ==%
        S = savingspline( mSavingsPar(:,ieps) );
        mAssetsPrime(:,ieps) = interp_savspline(S, MP.AssetsGridQuad);

    end

    % Initialize iteration
    err = 100; iter = 1;

    %== Options tested on solveStSt ==%
    options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
            'MaxFunEvals',50000,'TolX', 1e-16, 'TolFun', 1e-20, 'GradObj','on','MaxIter',5000);

    % options = optimoptions(@fminunc,'Algorithm','trust-region','Display','notify-detailed',...
    %         'MaxFunEvals',50000,'TolX', 1e-12,'TolFun', 1e-20, 'MaxIter',1000, ...
    %         'GradObj','on'); %,'DerivativeCheck','on','FinDiffType', 'central');

    mDensityCoeff = zeros(nMoments+1, neps);
    mDensity    = zeros(nAssetsQuad, neps);

    mMomentsCheck   = zeros(nMoments, neps);
    mMomentsNew     = zeros(nMoments, neps);
    aGridMomentsNew = zeros(nAssetsQuad, nMoments);

    while err > tolmoments && iter <= maxiter

        %== Find parameters of the density for given moments ==%
        for ieps = 1 : neps
            objectiveFunction = @(vParametersTilde) parametersResidual( vParametersTilde, aGridMoments(:,:,ieps) );

            [vParameters,normalization] = fminunc(objectiveFunction, zeros(nMoments,1), options);

            % storage parameters
            mDensityCoeff(:,ieps) = [1 / normalization; vParameters];

            % storage density values
            mDensity(:,ieps) = 1 / normalization * exp( aGridMoments(:,:,ieps) * vParameters );

            % CHECK moments × coefficients
            mMomentsCheck(1,ieps) = MP.QuadWeights' * ( MP.AssetsGridQuad .* mDensity(:,ieps) );

            for iMoment = 2 : nMoments
                mMomentsCheck(iMoment,ieps) = MP.QuadWeights' * ( ...
                        (MP.AssetsGridQuad - mMomentsCheck(1,ieps) ) .^ iMoment .* ...
                        mDensity(:,ieps) );
            end
        end

        %== CONSITENCY CHECK ==%
        err_par = max( abs( mMoments(:) - mMomentsCheck(:) ));

        try
            % if < 1e-6 it will cause trouble ahead when I check the STST
            assert( err_par < toldensity, 'parametersResidual couldnt`t find solution ')
        catch ME

            if print
                fprintf(1,'     Iteration %d Moment CHECK = %0.3e;  \n', iter, err_par );
            end
            % [mMoments(:) mMomentsCheck(:)];

            %options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','final',...
            %'MaxFunEvals',50000,'TolX', 1e-16, 'TolFun', 1e-20, 'GradObj','on','MaxIter',5000);

            options = optimoptions(@fminunc,'Algorithm','trust-region','Display','notify-detailed',...
            'MaxFunEvals',50000,'TolX', 1e-12,'TolFun', 1e-20, 'MaxIter',1000, ...
            'GradObj','on');

            for ieps = 1:neps
                objectiveFunction = @(vParametersTilde) parametersResidual( vParametersTilde, aGridMoments(:,:,ieps) );
                [vParametersNew,normalizationNew] = fminunc(objectiveFunction, vParameters, options);

                mDensityCoeff(:,ieps) = [1 / normalizationNew; vParametersNew];

                % storage density values
                mDensity(:,ieps) = 1 / normalizationNew * exp( aGridMoments(:,:,ieps) * vParametersNew );

                % CHECK moments Ã— coefficients
                mMomentsCheck(1,ieps) = MP.QuadWeights' * ( MP.AssetsGridQuad .* mDensity(:,ieps) );

               for iMoment = 2 : nMoments
                    mMomentsCheck(iMoment,ieps) = MP.QuadWeights' * ( ...
                           (MP.AssetsGridQuad - mMomentsCheck(1,ieps) ) .^ iMoment .* ...
                           mDensity(:,ieps) );
                end
            end
            [mMoments(:) mMomentsCheck(:)];
            err_parNew = max( abs( mMoments(:) - mMomentsCheck(:) ));

            if print
                fprintf(1,'     Iteration %d Moment Re-CHECK = %0.3e;  \n', iter, err_parNew );
            end
            assert( err_parNew < 10*toldensity, 'parametersResidual couldnt`t find solution AGAIN')
        end

        %== CONSITENCY CHECK ==%
        % err_par = max( abs( mMoments(:) - mMomentsCheck(:) ));
        % [mMoments(:) mMomentsCheck(:)]
        % fprintf(1,'     Moment CHECK = %0.5f;  \n', max( abs( mMoments(:) - mMomentsCheck(:) )) );

        %---------------------------------------------------------------%
        %%%  Compute next period moments and centered moments grid    %%%
        %---------------------------------------------------------------%
        for ieps = 1 : neps % WARN NEXT period shock

            %%%  NOTE:  computing next period first moments     %%%
            %%% m' = ∑ π( ̃ϵ | ϵ ) ∫ a'( ̃ϵ , a; z, m ) g_ϵ(a) da
            LawofMotionAux = 0;

            for jeps = 1 : neps % WARN THIS period shock

                LawofMotionAux = LawofMotionAux + ...
                        MP.ydist(jeps) * MP.transpp(jeps, ieps) * ...
                        MP.QuadWeights' * ( mAssetsPrime(:,jeps) .* mDensity(:, jeps) );
            end
            mMomentsNew(1,ieps) = LawofMotionAux / MP.ydist(ieps);

            aGridMomentsNew(:,1,ieps) = MP.AssetsGridQuad - mMomentsNew(1, ieps);

            % Compute higher order moments (centered)
            for iMoment = 2 : nMoments

                %%%  NOTE:  computing next period ith moment                    %%%
                %%%         Moment inside is with respect to ieps               %%%
                %%%
                %%%     mᵢ'(ϵ) = ∑ π( ̃ϵ | ϵ ) ∫ [ a'( ̃ϵ , a; z, m ) - m¹'(ϵ) ]^i g_̃ϵ (a) da
                %%%
                LawofMotionAux = 0;

                for jeps = 1 : neps % WARN THIS period shock

                    LawofMotionAux = LawofMotionAux + ...
                        MP.ydist(jeps) * MP.transpp(jeps, ieps) * ...
                        MP.QuadWeights' * ( ...
                        ( ( mAssetsPrime(:, jeps) - mMomentsNew(1,ieps) ) .^ iMoment ) .* ...
                        mDensity(:,jeps) );
                end

                mMomentsNew(iMoment, ieps) = LawofMotionAux / MP.ydist(ieps);

                aGridMomentsNew(:,iMoment, ieps) = ( MP.AssetsGridQuad - mMomentsNew(1,ieps) ) .^ iMoment - ...
                    mMomentsNew(iMoment, ieps);
            end

        end

        % **************************************************************************************

        %== Distance between moments ==%
        err = max( abs( mMomentsNew(:) - mMoments(:) ) );

        % print iterations.
        if mod(iter,50) ==0 && print
            fprintf('    Moments it %3d distance %.3e \n', iter, err)
            % fprintf('.\n   ')
        elseif print
            % fprintf('.')
        end
        iter = iter + 1;

        %== UPDATE ==%
        mMoments = mMomentsNew;
        aGridMoments = aGridMomentsNew;

        mMomentsNew(:)      = 0;
        aGridMomentsNew(:)  = 0;

        if iter == maxiter+1
            fprintf(1,'STOP at maxiter - MOMENT dist = %0.5f\n', err);

        end


    end

    % residual
    Ksupply = mMoments(1,:) * MP.ydist;
    resid = Ksupply - K;

    if print
        fprintf('\n\n')
        fprintf(1,'  K resid = %0.5f;  \n', resid);
        disp(['........................................................................'])
    end

    if nargout > 1
        mSavingsParOptional   = mSavingsPar;
        mDensityCoeffOptional = mDensityCoeff;
        mMomentsOptional      = mMoments;

        if nargout>4

            %== States ==%
            SS.mMoments         = mMoments;
            SS.mDensityCoeff    = mDensityCoeff;
            SS.mDensityExpCoeff = mDensityCoeff(2:end,:);

            SS.K             = Ksupply;

            %== policy ==%
            SS.mSavingsPar = mSavingsPar;

            %== Prices ==%
            SS.R = R;
            SS.wage = wage;
        end
    end
end


%%% Description:
%       Evaluates the euler residuals at stst
function res = eulerres_stst(vSavingsPar, R, wage)

    res = eulerres(vSavingsPar, vSavingsPar, R, R, wage, wage);
end
