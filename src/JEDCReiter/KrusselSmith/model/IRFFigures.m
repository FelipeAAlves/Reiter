


%== Histogram IRF Series ==%
vSavingsParSS = SS_Histogram.mSavingsPar(:);
irfHist_SavingsPar = irfHist(:, iSeriesHist.SavingsPar )' + repmat(vSavingsParSS, [1, MP.IRF_Length]);

irfHist_HistogramX = irfHist(:, iSeriesHist.HistogramX )';
irfHist_Histogram     = zeros(MP.nHistogramTotal, MP.IRF_Length);
irfHist_HistogramAux  = zeros(MP.nHistogramTotal, MP.IRF_Length);

irfHist_HistogramAux(:,1) = SS_Histogram.vHistogram;
for t = 1 : MP.IRF_Length
    %== Histogram from Solution ==%
    irfHist_Histogram(:,t) = x2distr(irfHist_HistogramX(:,t), SS_Histogram.vHistogram);

    %== Construct Histogram ==%
    if t>1
        %== Create Transition ==%
        msavingsPar = reshape(irfHist_SavingsPar(:,t-1), MP.nSavingsPar, MP.neps );
        Pi = sparse( forwardmat(0, msavingsPar) );

        %== Iterate foward ==%
        irfHist_HistogramAux(:,t) = forward(irfHist_HistogramAux(:,t-1),Pi);
    end
end

irfHist_Aggr = irfHist(:, iSeriesHist.Aggr );

% Color = [8/255*ones(MP.IRF_Length,1) 62/255*ones(MP.IRF_Length,1) (1-205/255*[1:MP.IRF_Length]'/MP.IRF_Length)];
% Color = cool( MP.IRF_Length );
% Color = spring(MP.IRF_Length);
Color = bone(MP.IRF_Length);
% Color = summer(MP.IRF_Length);
% Color = autumn(MP.IRF_Length);
% Color = winter(MP.IRF_Length);
% Color = lines(MP.IRF_Length);
% Color = gray(MP.IRF_Length);
% set(gca, 'ColorOrder', Color, 'NextPlot', 'replacechildren');

% ******************************************************************
%   FIGURES
% ==================================================================
%%
nHist = MP.nHistogram;
figure
AX(1) = subplot(1,2,1);
hold all
plot(MP.AssetsGridFine, irfHist_Histogram(1:nHist,1) / sum( irfHist_Histogram(1:nHist,1) ) ,'linewidth',3.0,...
'Color', [178/255, 34/255, 34/255],'DisplayName',sprintf('IRF %2d',1))
for t = MP.IRF_Length:-4:2
    plot(MP.AssetsGridFine, irfHist_Histogram(1:nHist,t) / sum( irfHist_Histogram(1:nHist,t) ) ,'linewidth',2.0,...
    'Color', Color(t,:))%,'DisplayName',sprintf('IRF %2d',t))
end
title('Unemployed Distribution IRF (Method)')
xlim([0,9])
hold off
ylim1 = ylim;

subplot(1,2,2);
hold all
plot(MP.AssetsGridFine, irfHist_Histogram(1:nHist,1) / sum( irfHist_Histogram(1:nHist,1) ) ,'linewidth',3.0,...
'Color', [178/255, 34/255, 34/255],'DisplayName',sprintf('IRF %2d',1))
for t = MP.IRF_Length:-4:2
    plot(MP.AssetsGridFine, irfHist_HistogramAux(1:nHist,t) / sum( irfHist_HistogramAux(1:nHist,t) ) ,'linewidth',2.0,...
        'color', Color(t,:), 'DisplayName',sprintf('IRF %2d',t))
end
legend('show')
title('Unemployed Distribution IRF (Auxiliary)')
xlim([0,9])
ylim(ylim1)
hold off

allYLim = get(AX, {'YLim'});
allYLim = cat(2, allYLim{:});
set(AX, 'YLim', [min(allYLim), max(allYLim)]);
% --------------------------------------------------------------------------------------
% % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
% --------------------------------------------------------------------------------------

figure
subplot(1,2,1)
hold all
plot(MP.AssetsGridFine, irfHist_Histogram(nHist+1:end,1) / sum( irfHist_Histogram(nHist+1:end,1) ) ,'linewidth',3.0,...
'Color', [178/255, 34/255, 34/255],'DisplayName',sprintf('IRF %2d',1))
% plot(MP.AssetsGridFine, irfHist_Histogram(nHist+1:end,2) / sum( irfHist_Histogram(nHist+1:end,2) ) ,'linewidth',3.0,...
% 'Color', [125/255, 34/255, 34/255],'DisplayName',sprintf('IRF %2d',2))
% plot(MP.AssetsGridFine, irfHist_Histogram(nHist+1:end,3) / sum( irfHist_Histogram(nHist+1:end,3) ) ,'linewidth',3.0,...
% 'Color', [75/255, 34/255, 34/255],'DisplayName',sprintf('IRF %2d',3))
% plot(MP.AssetsGridFine, irfHist_Histogram(nHist+1:end,4) / sum( irfHist_Histogram(nHist+1:end,4) ) ,'linewidth',3.0,...
% 'Color', [25/255, 34/255, 34/255],'DisplayName',sprintf('IRF %2d',4))
for t = MP.IRF_Length:-4:2
    plot(MP.AssetsGridFine, irfHist_Histogram(nHist+1:end,t) / sum( irfHist_Histogram(nHist+1:end,t) ) ,'linewidth',2.0,...
    'Color', Color(t,:),'DisplayName',sprintf('IRF %2d',t))
end
legend('show')
title('Employed Distribution IRF (Method)')
xlabel('Assets, $a$','interpreter','latex')
ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
xlim([MP.aabar .5*MP.AssetsMax])
set(gcf,'color','w')
hold off
ylim1 = ylim;

% .....................................................................................
subplot(1,2,2)
hold all
plot(MP.AssetsGridFine, irfHist_Histogram(nHist+1:end,1) / sum( irfHist_Histogram(nHist+1:end,1) ) ,'linewidth',3.0,...
'Color', [178/255, 34/255, 34/255],'DisplayName',sprintf('IRF %2d',1))
for t = MP.IRF_Length:-4:2
    plot(MP.AssetsGridFine, irfHist_HistogramAux(nHist+1:end,t) / sum( irfHist_HistogramAux(nHist+1:end,t) ) ,'linewidth',2.0,...
        'color', Color(t,:), 'DisplayName',sprintf('IRF %2d',t))
end
legend('show')
title('Employed Distribution IRF (Auxiliary)')
xlabel('Assets, $a$','interpreter','latex')
ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
ylim(ylim1)
xlim([MP.aabar .5*MP.AssetsMax])
set(gcf,'color','w')
hold off
