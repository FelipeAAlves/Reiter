% ******************************************************************
%   COMPUTE Steady-State
% ==================================================================

fprintf('Compute steady state from Histogram...\n')

Kstart = [MP.KRepSS; 2*MP.KRepSS];
% Kstart = [1.0; 1.8];

%%% WINBERRY fsolve options
% options = optimoptions('fsolve','Display','iter-detailed');
%%% REITER fzero options
opts = optimset('Display', 'off', 'TolX',  1e-10);

%== Find Kstst ==%
KequHistogram = fzero(@ststHistogramResid, Kstart, opts);

%== Recover everything ==%
[resid, mSavingsPar, mHistogram, SS_Histogram] = ststHistogramResid(KequHistogram);
fprintf(1,'\nStSt done: Kequ = %0.6f\n',KequHistogram);

%%% RESULTS
%%% np=3 - 1.620279

% ******************************************************************
%   FIGURES
% ==================================================================
if plt_figure

    mSavingsHist = zeros(MP.nHistogram, MP.neps);
    mConsumptionHist = zeros(MP.nHistogram, MP.neps);
    % mSavingsPar     = SS_Histogram.mSavingsPar;

    for ieps = 1 : MP.neps

        S = savingspline( mSavingsPar(:,ieps) );
        mSavingsHist(:,ieps) = interp_savspline(S,MP.AssetsGridFine);

        %== Prices ==%
        Rstst      = SS_Histogram.R;
        wagestst   = SS_Histogram.wage;

        mConsumptionHist(:, ieps) = getConsumption( ...
                    MP.AssetsGridFine, mSavingsHist(:, ieps), Rstst, wagestst, ieps);

    end

    if Winberry == 0
        figure
        hold on
        plot(MP.AssetsGridFine, mHistogram(:,1) / sum( mHistogram(:,1) ),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mHistogram(:,2) / sum( mHistogram(:,2) ),'linewidth',1.5,'color',[178/255,34/255,34/255])
        xlabel('Assets, $a$', 'interpreter', 'latex')
        ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
        xlim([MP.aabar .9*MP.AssetsMax])
        title('Invariant Distribution of Households')
        legend('Histogram Unemployed','Histogram Employed','location','northeast')
        % set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
        set(gcf,'color','w')
        hold off
    end
end
%======================================%
%%% CHECK Euler resid on Stst        %%%
%======================================%
% Par(1) is kink point of consumption function_:
% xx2 = linspace(Par(1),MP.xmax,10000)';
% R = 1+netintr(Kequ,1)-MP.taustst;
% wage = wagefunc(Kequ,1);
% Transf = MP.taustst*Kequ;
% resC = eulerres(Par,Par,R,wage,Transf,xx2);
% MP.eulerres = [max(abs(resC)) mean(abs(resC))];

% ******************************************************************
%   Parametric density StSt
% ==================================================================
if Winberry > 0

    fprintf('\n')
    fprintf('Compute steady state from parametric family...\n\n')

    % Compute moments from histogram
    mMomentsHistogram  = zeros(MP.nMoments, MP.neps);
    aGridMoments       = zeros(MP.nAssetsQuad, MP.nMoments, MP.neps);

    % CHECK the parametersResidual fnc
    options1 = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
            'MaxFunEvals',50000,'TolX', 1e-12, 'TolFun', 1e-15, 'GradObj','on','MaxIter',1000);

    options2 = optimoptions(@fminunc,'Algorithm','trust-region','Display','notify-detailed',...
            'MaxFunEvals',50000,'TolX', 1e-12,'TolFun', 1e-15, 'MaxIter',1000, ...
            'GradObj','on'); %,'DerivativeCheck','on','FinDiffType', 'central');

    mDensity1       = zeros(MP.nAssetsQuad, MP.neps);
    mMomentsCheck1  = zeros(MP.nMoments, MP.neps);

    mDensity2        = zeros(MP.nAssetsQuad, MP.neps);
    mMomentsCheck2   = zeros(MP.nMoments, MP.neps);

    aFineGridMoments   = zeros(MP.nHistogram, MP.nMoments, MP.neps);
    mDensityFine    = zeros(MP.nHistogram, MP.neps);

    for ieps = 1 : MP.neps

        % First moment (uncentered)
        mMomentsHistogram(1,ieps) = dot( MP.AssetsGridFine , ( mHistogram(:, ieps) ./ sum(mHistogram(:,ieps)) ) );

        aGridMoments(:,1, ieps)     = MP.AssetsGridQuad - mMomentsHistogram(1,ieps);
        aFineGridMoments(:,1, ieps) = MP.AssetsGridFine - mMomentsHistogram(1,ieps);

        % Higher order moments (centered)
        for iMoment = 2 : MP.nMoments
            mMomentsHistogram(iMoment,ieps) = dot(...
                ( MP.AssetsGridFine - mMomentsHistogram(1,ieps) ) .^ iMoment, ...
                mHistogram(:, ieps) ./ sum( mHistogram(:,ieps) ) ...
            );

            % aMoments on the QuadNodes
            aGridMoments(:,iMoment,ieps) = ...
            (MP.AssetsGridQuad - mMomentsHistogram(1,ieps) ) .^ iMoment - mMomentsHistogram(iMoment,ieps);

            % aMoments on the FineGrid
            aFineGridMoments(:,iMoment,ieps) = ...
            (MP.AssetsGridFine - mMomentsHistogram(1,ieps) ) .^ iMoment - mMomentsHistogram(iMoment,ieps);

        end

        objectiveFunction = @(vParametersTilde) parametersResidual( vParametersTilde, aGridMoments(:,:,ieps) );


        %%%  NOTE:  Compare how are two OPTIONS %%%
        % g = @() fminunc(objectiveFunction, zeros(MP.nMoments,1), options1);
        % h = @() fminunc(objectiveFunction, zeros(MP.nMoments,1), options2);
        % fprintf('option0 : %.6f , option : %.6f\n', timeit(g), timeit(h))

        [vParameters1,normalization1] = fminunc(objectiveFunction, zeros(MP.nMoments,1), options1);
        [vParameters2,normalization2] = fminunc(objectiveFunction, zeros(MP.nMoments,1), options2);

        % storage density values
        mDensity1(:,ieps) = 1 / normalization1 * exp( aGridMoments(:,:,ieps) * vParameters1 );
        mDensity2(:,ieps)  = 1 / normalization2 * exp( aGridMoments(:,:,ieps) * vParameters2 );

        % storage density values on FineGrid
        mDensityFine(:,ieps) = 1 / normalization1 * exp( aFineGridMoments(:,:,ieps) * vParameters1 );

        % CHECK moments ?? coefficients
        mMomentsCheck1(1,ieps) = MP.QuadWeights' * ( MP.AssetsGridQuad .* mDensity1(:,ieps) );
        mMomentsCheck2(1,ieps) = MP.QuadWeights' * ( MP.AssetsGridQuad .* mDensity2(:,ieps) );

        for iMoment = 2 : MP.nMoments
            mMomentsCheck1(iMoment,ieps) = MP.QuadWeights' * ( ...
            (MP.AssetsGridQuad - mMomentsCheck1(1,ieps) ) .^ iMoment .* ...
            mDensity1(:,ieps) );

            mMomentsCheck2(iMoment,ieps) = MP.QuadWeights' * ( ...
            (MP.AssetsGridQuad - mMomentsCheck2(1,ieps) ) .^ iMoment .* ...
            mDensity2(:,ieps) );
        end

    end

    err_par1 = max( abs( mMomentsHistogram(:) - mMomentsCheck1(:) ));
    err_par2 = max( abs( mMomentsHistogram(:) - mMomentsCheck2(:) ));

    fprintf('Moments Histogram x Moments Density \n')
    disp([mMomentsHistogram(:) mMomentsCheck1(:) mMomentsCheck2(:)])
    fprintf('distance for quasi-newton algorthim: %.4e\n', err_par1)
    fprintf('distance for trust-region algorthim: %.4e\n', err_par2)
    disp(['........................................................................'])

    % Histogram × parametric
    if plt_figure
        figure
        hold on
        plot(MP.AssetsGridFine, mHistogram(:,1) / sum( mHistogram(:,1) )  ,'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mDensityFine(:,1) / sum(mDensityFine(:,1)),'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
        xlabel('Assets, $a$','interpreter','latex')
        ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
        xlim([MP.aabar .9*MP.AssetsMax])
        title('Invariant Distribution of Households (Unemployed)')
        legend('Histogram','Parametric Family','location','northeast')
        set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
        set(gcf,'color','w')
        hold off
    end
    % .....................................................................................

    %== Evaluate residPolynom at previous equilibrium ==%
    fprintf('EVAL residPolynomial at previous Kequ\n')
    residPolynom = ststPolynomialResid(KequHistogram, mMomentsHistogram, aGridMoments,1);
    fprintf('ststPolynomialResid at previous Kequ = %.6f\n', residPolynom )

    if abs(residPolynom) > 1e-4

        f = @(K) ststPolynomialResid(K,mMomentsHistogram,aGridMoments);
        % opts = optimset('Display', 'off', 'TolX',  1e-10);
        % KequPolynomial = fzero(f, KequHistogram, opts);

        %%%  NOTE:  fsolve performs better at this stage     %%%
        options = optimoptions('fsolve','Display','iter-detailed');
        [KequPolynomial,err,exitflag] = fsolve(f,KequHistogram,options);
    end

    [resid, mSavingsPar, mDensityCoeff, mMoments, SS_Polynomial] = ...
    ststPolynomialResid(KequPolynomial, mMomentsHistogram, aGridMoments);
    % ******************************************************************
    %   FIGURES
    % ==================================================================
    if plt_figure

        %% Compute exponential polynomial pdf along fine grid %%
        % *********************************************************************

        mDensityFine = zeros(MP.nHistogram, MP.neps);
        aGridFineMoments = zeros(MP.nHistogram, MP.nMoments);

        for ieps = 1 : MP.neps

            % First moment (uncentered) - analogous to aGridMoments
            aGridFineMoments(:,1) = ( MP.AssetsGridFine - mMoments(1,ieps) );

            % Higher order moments (centered)
            for iMoment = 2 : MP.nMoments
                aGridFineMoments(:,iMoment) = (MP.AssetsGridFine - mMoments(1,ieps)) .^ iMoment - ...
                mMoments(iMoment,ieps);
            end

            % Compute distribution
            mDensityFine(:, ieps) = mDensityCoeff(1,ieps) * exp(aGridFineMoments * ...
            mDensityCoeff(2:end, ieps) );

        end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Compute policy over fine Grid %%
        % *********************************************************************

        mSavingsPoly = zeros(MP.nHistogram, MP.neps);
        mConsumptionPoly = zeros(MP.nHistogram, MP.neps);
        % mSavingsPar     = SS_Polynomial.mSavingsPar;

        for ieps = 1 : MP.neps

            %== savings policy  ==%
            S = savingspline( mSavingsPar(:,ieps) );

            mSavingsPoly(:, ieps) = interp_savspline(S,MP.AssetsGridFine);

            %== Prices ==%
            Rstst      = SS_Polynomial.R;
            wagestst   = SS_Polynomial.wage;

            mConsumptionPoly(:, ieps) = getConsumption( ...
                        MP.AssetsGridFine, mSavingsPoly(:, ieps), Rstst, wagestst, ieps);
        end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Policies
        figure
        hold on
        plot(MP.AssetsGridFine, mSavingsHist(:,1),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mSavingsPoly(:,1),'linewidth',1.5,'color',[178/255,34/255,34/255])
        plot(MP.AssetsGridFine, mSavingsHist(:,2),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mSavingsPoly(:,2),'linewidth',1.5,'color',[178/255,34/255,34/255])
        xlabel('Assets, $a$','interpreter','latex')
        ylabel('Savings, $a`(a, \varepsilon)$','interpreter','latex')
        xlim([MP.aabar .25*MP.AssetsMax])
        title('Savings Policy')
        legend('Hist Unemp','Poly Unemp','Hist Empl','Poly Empl','location','northeast')
        % set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
        set(gcf,'color','w')
        hold off

        figure
        subplot(1,2,1)
        hold on
        plot(MP.AssetsGridFine, mConsumptionHist(:,1),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mConsumptionPoly(:,1),'linewidth',1.5,'color',[178/255,34/255,34/255])
        plot(MP.AssetsGridFine, mConsumptionHist(:,2),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mConsumptionPoly(:,2),'linewidth',1.5,'color',[178/255,34/255,34/255])
        xlabel('Assets, $a$','interpreter','latex')
        ylabel('Consumption, $c(a, \varepsilon)$','interpreter','latex')
        xlim([MP.aabar .25*MP.AssetsMax])
        title('Consumption Policy')
        legend('Hist Unemp','Poly Unemp','Hist Empl','Poly Empl','location','northeast')
        % set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
        set(gcf,'color','w')
        hold off
        % .....................................................................................
        subplot(1,2,2)
        hold on
        plot(MP.AssetsGridFine, mConsumptionHist(:,1),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mConsumptionPoly(:,1),'linewidth',1.5,'color',[178/255,34/255,34/255])
        plot(MP.AssetsGridFine, mConsumptionHist(:,2),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mConsumptionPoly(:,2),'linewidth',1.5,'color',[178/255,34/255,34/255])
        xlabel('Assets, $a$','interpreter','latex')
        ylabel('Consumption, $c(a, \varepsilon)$','interpreter','latex')
        xlim([MP.aabar .8*MP.AssetsMax])
        title('Consumption Policy')
        legend('Hist Unemp','Poly Unemp','Hist Empl','Poly Empl','location','northeast')
        % set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
        set(gcf,'color','w')
        hold off

        % --------------------------------------------------------------------------------------
        % % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
        % --------------------------------------------------------------------------------------

        % Histogram ?? Density
        figure
        hold on
        plot(MP.AssetsGridFine, mHistogram(:,1) / sum( mHistogram(:,1) ),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mDensityFine(:,1) / sum(mDensityFine(:,1)),'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
        xlabel('Assets, $a$','interpreter','latex')
        ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
        xlim([MP.aabar 1.0*MP.AssetsMax])
        title('Invariant Distribution of Households (Unemployed)')
        legend('Histogram','Parametric Family','location','northeast')
        set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
        set(gcf,'color','w')
        hold off

        figure
        hold on
        plot(MP.AssetsGridFine, mHistogram(:,2) / sum( mHistogram(:,2) ),'linewidth',1.5,'color',[8/255,62/255,118/255])
        plot(MP.AssetsGridFine, mDensityFine(:,2) / sum(mDensityFine(:,2)),'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
        xlabel('Assets, $a$','interpreter','latex')
        ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
        xlim([MP.aabar 1.0*MP.AssetsMax])
        title('Invariant Distribution of Households (Employed)')
        legend('Histogram','Parametric Family','location','northeast')
        set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
        set(gcf,'color','w')
        hold off
    end
end
