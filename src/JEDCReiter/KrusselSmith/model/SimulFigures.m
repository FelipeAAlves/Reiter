%% Work with SIMULATED Series %%
% *********************************************************************

if Reiter
    %== Histogram Sim Series ==%
    simulHist_HistogramX = simulHist(:, iSeriesHist.HistogramX )';
    simulHist_Histogram  = zeros(MP.nHistogramTotal, MP.sim_T);
    for t = 1 : MP.sim_T
        simulHist_Histogram(:,t) = x2distr(simulHist_HistogramX(:,t), SS_Histogram.vHistogram);
    end

    simulHist_Aggr = simulHist(:, iSeriesHist.Aggr );

    vSavingsParSS = SS_Histogram.mSavingsPar(:);
    simulHist_SavingsPar = simulHist(:, iSeriesHist.SavingsPar )' + repmat(vSavingsParSS, [1, MP.sim_T]);

    ZpltHist = simulHist_Aggr(:,1);
    KpltHistogram  = simulHist_Aggr(:,2) + SS_Histogram.K;

    %--------------------------------%
    %%% High and LOW period        %%%
    %--------------------------------%
    sig = sqrt(MP.sigTFP^2/(1-MP.rhoTFP^2));
    [~,indHigh] = max(ZpltHist);
    [~,indLow ] = min(ZpltHist);
    indHigh = indHigh + 1;
    indLow  = indLow  + 1;
%     indHigh = find(ZpltHist >  3*sig,1) + 1;    % NEXT period to High shock
%     indLow  = find(ZpltHist < -3*sig,1) + 1;

    % histogram
    mHistogramSS   = reshape(SS_Histogram.vHistogram, MP.nHistogram, MP.neps);
    mHistogramHigh = reshape(simulHist_Histogram(:, indHigh), MP.nHistogram, MP.neps);
    mHistogramLow  = reshape(simulHist_Histogram(:, indLow) , MP.nHistogram, MP.neps);

    %------------------------------------------%
    %%% Policies - for Histogram only        %%%
    %------------------------------------------%
    mSavingsParSS   = SS_Histogram.mSavingsPar;
    mSavingsParHigh = reshape( simulHist_SavingsPar(:,indHigh-1), MP.nSavingsPar, MP.neps);
    mSavingsParLow  = reshape( simulHist_SavingsPar(:,indLow-1) , MP.nSavingsPar, MP.neps);

    mSavingsSS   = zeros(MP.nHistogram, MP.neps);
    mSavingsHigh = zeros(MP.nHistogram, MP.neps);
    mSavingsLow  = zeros(MP.nHistogram, MP.neps);


    mConsumptionSS   = zeros(MP.nHistogram, MP.neps);
    mConsumptionHigh = zeros(MP.nHistogram, MP.neps);
    mConsumptionLow  = zeros(MP.nHistogram, MP.neps);

    for ieps = 1 : MP.neps


        %%%     Steady State     %%%
        S = savingspline( mSavingsParSS(:,ieps) );
        mSavingsSS(:,ieps) = interp_savspline(S,MP.AssetsGridFine);

        %== Prices ==%
        Rstst      = SS_Histogram.R;
        wagestst   = SS_Histogram.wage;

        mConsumptionSS(:, ieps) = getConsumption( ...
                    MP.AssetsGridFine, mSavingsSS(:, ieps), Rstst, wagestst, ieps);

        % ....................................................................................

        %%%     HIGH     %%%
        S = savingspline( mSavingsParHigh(:,ieps) );
        mSavingsHigh(:,ieps) = interp_savspline(S,MP.AssetsGridFine);

        K_High = simulHist_Aggr(indHigh-1,2) + SS_Histogram.K;
        dzHigh = simulHist_Aggr(indHigh-1,1) ;

        RHigh      = 1 + netintr(K_High, 1 + dzHigh);
        wageHigh   = wagefunc(K_High, 1 + dzHigh);

        mConsumptionHigh(:, ieps) = getConsumption( ...
        MP.AssetsGridFine, mSavingsHigh(:, ieps), RHigh, wageHigh, ieps);

        % ....................................................................................

        %%%     LOW     %%%
        S = savingspline( mSavingsParLow(:,ieps) );
        mSavingsLow(:,ieps) = interp_savspline(S,MP.AssetsGridFine);

        K_Low = simulHist_Aggr(indLow-1,2) + SS_Histogram.K;
        dzLow = simulHist_Aggr(indLow-1,1) ;

        RLow      = 1 + netintr(K_Low, 1 + dzLow);
        wageLow   = wagefunc(K_Low, 1 + dzLow);

        mConsumptionLow(:, ieps) = getConsumption( ...
                    MP.AssetsGridFine, mSavingsLow(:, ieps), RLow, wageLow, ieps);

    end
end
% .....................................................................................
if Winberry

    %== Polynomial Series ==%
    vMomentsSS = SS_Polynomial.mMoments(:);
    vDensityExpCoeffSS = SS_Polynomial.mDensityExpCoeff(:);

    simPolyMoments         = seriesPolynomial(:,iSeriesPoly.Mom)'  + repmat( vMomentsSS, [1, MP.sim_T] );
    simPolyDensityExpCoeff = seriesPolynomial(:,iSeriesPoly.Dens)' + repmat( vDensityExpCoeffSS, [1, MP.sim_T] );

    simPolyDensityCoeff  = zeros(MP.nDensityCoeffTotal, MP.sim_T);
    for t = 1 : MP.sim_T
        [~, mDensityCoeff] = XtoDensity( simPolyMoments(:,t), simPolyDensityExpCoeff(:,t) );
        simPolyDensityCoeff(:,t) = mDensityCoeff(:);
    end

    simPolyAggr = seriesPolynomial(:,iSeriesPoly.Aggr);


    ZpltPoly = simPolyAggr(:,1);

    KpltPolynomial = simPolyAggr(:,2) + SS_Polynomial.K;

    %-------------------------%
    %%% Distribution        %%%
    %-------------------------%

    % density coeff
    mDensityCoeffSS   = SS_Polynomial.mDensityCoeff;
    mDensityCoeffHigh = reshape(simPolyDensityCoeff(:,indHigh), MP.nDensityCoeff, MP.neps);
    mDensityCoeffLow  = reshape(simPolyDensityCoeff(:,indLow) , MP.nDensityCoeff, MP.neps);

    % moments
    mMomentsSS   = SS_Polynomial.mMoments;
    mMomentsHigh = reshape( simPolyMoments(:, indHigh), MP.nMoments, MP.neps);
    mMomentsLow  = reshape( simPolyMoments(:, indLow) , MP.nMoments, MP.neps);

    %== Pre-Allocate ==%
    aFineGridMomentsSS   = zeros(MP.nHistogram, MP.nMoments, MP.neps);
    aFineGridMomentsHigh = zeros(MP.nHistogram, MP.nMoments, MP.neps);
    aFineGridMomentsLow  = zeros(MP.nHistogram, MP.nMoments, MP.neps);

    mDensityFineSS   = zeros(MP.nHistogram, MP.neps);
    mDensityFineHigh = zeros(MP.nHistogram, MP.neps);
    mDensityFineLow  = zeros(MP.nHistogram, MP.neps);

    for ieps =1 :MP.neps

        aFineGridMomentsSS(:,1, ieps)   = MP.AssetsGridFine - mMomentsSS(1,ieps);
        aFineGridMomentsHigh(:,1, ieps) = MP.AssetsGridFine - mMomentsHigh(1,ieps);
        aFineGridMomentsLow(:,1, ieps)  = MP.AssetsGridFine - mMomentsLow(1,ieps);

        % Higher order moments (centered)
        for iMoment = 2 : MP.nMoments


            aFineGridMomentsSS(:,iMoment,ieps) = ...
            ( MP.AssetsGridFine - mMomentsSS(1,ieps) ) .^ iMoment - mMomentsSS(iMoment,ieps);

            aFineGridMomentsHigh(:,iMoment,ieps) = ...
            ( MP.AssetsGridFine - mMomentsHigh(1,ieps) ) .^ iMoment - mMomentsHigh(iMoment,ieps);

            aFineGridMomentsLow(:,iMoment,ieps) = ...
            (MP.AssetsGridFine - mMomentsLow(1,ieps) ) .^ iMoment - mMomentsLow(iMoment,ieps);

        end

        % storage density values on FineGrid
        mDensityFineSS(:,ieps)   = 1 / mDensityCoeffSS(1,ieps) * exp( aFineGridMomentsSS(:,:,ieps) * mDensityCoeffSS(2:end,ieps) );
        mDensityFineHigh(:,ieps) = 1 / mDensityCoeffHigh(1,ieps) * exp( aFineGridMomentsHigh(:,:,ieps) * mDensityCoeffHigh(2:end,ieps) );
        mDensityFineLow(:,ieps)  = 1 / mDensityCoeffLow(1,ieps) * exp( aFineGridMomentsLow(:,:,ieps) * mDensityCoeffLow(2:end,ieps) );
    end

else
    ZpltPoly = zeros(100,1);
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ******************************************************************
%   FIGURES
% ==================================================================


%fprintf('Max difference between Histogram x Polynomial: %.5f\n', max(abs( KpltHistogram -  KpltPolynomial)) )

% Technological Shock
figure
hold on
plot(ZpltHist(1:100) , 'linewidth',2.5,'color',[178/255,34/255,34/255])
%plot(ZpltPoly(1:100) , 'linewidth',2.5,'color',[8/255,62/255,118/255])
xlabel('Time, $t$','interpreter','latex')
ylabel('Technology $z$','interpreter','latex')
% title('$K$ Simulation','interpreter','latex')
% set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
set(gcf,'color','w')
hold off

% Capital Figure
figure
hold on
plot(KpltHistogram(1:100) , 'linewidth',2.5,'color',[178/255,34/255,34/255])
%plot(KpltPolynomial(1:100), 'linewidth',2.5,'color',[8/255,62/255,118/255])
xlabel('Time, $t$','interpreter','latex')
ylabel('Capital')
title('$K$ Simulation','interpreter','latex')
legend('Histogram','Parametric Family','location','northeast')
% set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
set(gcf,'color','w')
hold off

% --------------------------------------------------------------------------------------
% % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
% --------------------------------------------------------------------------------------

figure
hold on
plot(MP.AssetsGridFine, mSavingsSS  (:,1),'linewidth',1.5,'color',[105/255,105/255,105/255])
plot(MP.AssetsGridFine, mSavingsSS  (:,2),'linewidth',1.5,'color',[105/255,105/255,105/255])
plot(MP.AssetsGridFine, mSavingsHigh(:,1),'linewidth',1.5,'color',[  8/255, 62/255,118/255])
plot(MP.AssetsGridFine, mSavingsHigh(:,2),'linewidth',1.5,'color',[  8/255, 62/255,118/255])
plot(MP.AssetsGridFine, mSavingsLow (:,1),'linewidth',1.5,'color',[178/255, 34/255, 34/255])
plot(MP.AssetsGridFine, mSavingsLow (:,2),'linewidth',1.5,'color',[178/255, 34/255, 34/255])
xlabel('Assets, $a$','interpreter','latex')
ylabel('Savings, $a`(a, \varepsilon)$','interpreter','latex')
xlim([MP.aabar .25*MP.AssetsMax])
title('Savings Policy')
% legend('Hist Unemp','Poly Unemp','Hist Empl','Poly Empl','location','northeast')
% set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
set(gcf,'color','w')
hold off

figure
subplot(1,2,1)
hold on
plot(MP.AssetsGridFine, mConsumptionSS  (:,1),'linewidth',1.5,'color',[105/255,105/255,105/255])
plot(MP.AssetsGridFine, mConsumptionSS  (:,2),'linewidth',1.5,'color',[105/255,105/255,105/255])
plot(MP.AssetsGridFine, mConsumptionHigh(:,1),'linewidth',1.5,'color',[  8/255, 62/255,118/255])
plot(MP.AssetsGridFine, mConsumptionHigh(:,2),'linewidth',1.5,'color',[  8/255, 62/255,118/255])
plot(MP.AssetsGridFine, mConsumptionLow (:,1),'linewidth',1.5,'color',[178/255, 34/255, 34/255])
plot(MP.AssetsGridFine, mConsumptionLow (:,2),'linewidth',1.5,'color',[178/255, 34/255, 34/255])
xlabel('Assets, $a$','interpreter','latex')
ylabel('Savings, $a`(a, \varepsilon)$','interpreter','latex')
xlim([MP.aabar .25*MP.AssetsMax])
title('Savings Policy')
% legend('Hist Unemp','Poly Unemp','Hist Empl','Poly Empl','location','northeast')
% set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
set(gcf,'color','w')
hold off
% .....................................................................................
subplot(1,2,2)
hold on
plot(MP.AssetsGridFine, mConsumptionSS  (:,1),'linewidth',1.5,'color',[105/255,105/255,105/255])
plot(MP.AssetsGridFine, mConsumptionSS  (:,2),'linewidth',1.5,'color',[105/255,105/255,105/255])
plot(MP.AssetsGridFine, mConsumptionHigh(:,1),'linewidth',1.5,'color',[  8/255, 62/255,118/255])
plot(MP.AssetsGridFine, mConsumptionHigh(:,2),'linewidth',1.5,'color',[  8/255, 62/255,118/255])
plot(MP.AssetsGridFine, mConsumptionLow (:,1),'linewidth',1.5,'color',[178/255, 34/255, 34/255])
plot(MP.AssetsGridFine, mConsumptionLow (:,2),'linewidth',1.5,'color',[178/255, 34/255, 34/255])
xlabel('Assets, $a$','interpreter','latex')
ylabel('Savings, $a`(a, \varepsilon)$','interpreter','latex')
xlim([MP.aabar .75*MP.AssetsMax])
title('Savings Policy')
% legend('Hist Unemp','Poly Unemp','Hist Empl','Poly Empl','location','northeast')
% set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
set(gcf,'color','w')
hold off

% --------------------------------------------------------------------------------------
% % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
% --------------------------------------------------------------------------------------

% Distribution High x Low shock Unemployed Figure
figure
hold on
plot(MP.AssetsGridFine, mHistogramSS(:,1) / sum( mHistogramSS(:,1) )      ,'linewidth',1.5,'color',[105/255,105/255,105/255])
%plot(MP.AssetsGridFine, mDensityFineSS(:,1) / sum(mDensityFineSS(:,1))    ,'linewidth',1.5,'color',[105/255,105/255,105/255],'linestyle','--')
plot(MP.AssetsGridFine, mHistogramHigh(:,1) / sum( mHistogramHigh(:,1) )  ,'linewidth',1.5,'color',[8/255  , 62/255,118/255])
%plot(MP.AssetsGridFine, mDensityFineHigh(:,1) / sum(mDensityFineHigh(:,1)),'linewidth',1.5,'color',[8/255  , 62/255,118/255],'linestyle','--')
plot(MP.AssetsGridFine, mHistogramLow(:,1) / sum( mHistogramLow(:,1) )    ,'linewidth',1.5,'color',[178/255, 34/255, 34/255])
%plot(MP.AssetsGridFine, mDensityFineLow(:,1) / sum(mDensityFineLow(:,1))  ,'linewidth',1.5,'color',[178/255, 34/255, 34/255],'linestyle','--')
xlabel('Assets, $a$','interpreter','latex')
ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
xlim([MP.aabar .9*MP.AssetsMax])
title('Distribution of Households (Unemployed)')
legend('Hist StSt','Dens StSt','Hist High','Dens High','Hist Low','Dens Low','location','best')
% set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
set(gcf,'color','w')
hold off

figure
hold on
plot(MP.AssetsGridFine, mHistogramSS(:,2) / sum( mHistogramSS(:,2) )      ,'linewidth',1.5,'color',[105/255,105/255,105/255])
%plot(MP.AssetsGridFine, mDensityFineSS(:,2) / sum(mDensityFineSS(:,2))    ,'linewidth',1.5,'color',[105/255,105/255,105/255],'linestyle','--')
plot(MP.AssetsGridFine, mHistogramHigh(:,2) / sum( mHistogramHigh(:,2) )  ,'linewidth',1.5,'color',[8/255  , 62/255,118/255])
%plot(MP.AssetsGridFine, mDensityFineHigh(:,2) / sum(mDensityFineHigh(:,2)),'linewidth',1.5,'color',[8/255  , 62/255,118/255],'linestyle','--')
plot(MP.AssetsGridFine, mHistogramLow(:,2) / sum( mHistogramLow(:,2) )    ,'linewidth',1.5,'color',[178/255, 34/255,34/255 ])
%plot(MP.AssetsGridFine, mDensityFineLow(:,2) / sum(mDensityFineLow(:,2))  ,'linewidth',1.5,'color',[178/255, 34/255,34/255 ],'linestyle','--')
xlabel('Assets, $a$','interpreter','latex')
ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
xlim([MP.aabar .9*MP.AssetsMax])
title('Invariant Distribution of Households (Employed)')
legend('Hist StSt','Dens StSt','Hist High','Dens High','Hist Low','Dens Low','location','best')
% set(gca,'xtick',[MP.aabar],'XTickLabel','$\underline{a}$','TickLabelInterpreter','latex')
set(gcf,'color','w')
hold off
