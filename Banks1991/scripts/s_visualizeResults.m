% Visualize results


% Load data

resultsDir = fullfile(banksRootPath, 'results');

whichObserver = 'ideal';
allSegments = {'inner', 'outer'};
thisEccen = 10;

for s = 1:length(allSegments)
    r(s) =  load(fullfile(saveDir, sprintf('dPrime_%s_%s_eccen%d.mat', whichObserver, allSegments{s}, thisEccen)));
    
end


% Load experiment parameters
expParams = loadExpParamsBanks1991;



%% debug ideal observer
threshold = 1.36/sqrt(2);

for ii = 1:9, label{ii} = sprintf('%1.3f',expParams.sf(1,ii)); end


cmap = jet(9);
for s = [1,2]
    figure; set(gcf, 'Color', 'w'); clf; hold on;
    for ii = 1:9
        plot(expParams.contrastLevels*100, r(s).thisdPrime(:,ii), 'Color', cmap(ii,:,:))
    end
    
    plot([0.1 100], [threshold, threshold],'k--')
    % plot([0.1 100], [1.36, 1.36],'r--')
    % xlim([0.1 100]); ylim([0.1 100]);
    xlabel('Contrast (Michelson)')
    ylabel('Contrast Sensitivity (d prime)'); title('Figure 3')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    set(gca, 'XTick', [0.1, 1, 10, 100], 'XTickLabel', {'0.1', '1', '10', '100'})
    set(gca, 'YTick', [0.1, 1, 10, 100], 'YTickLabel', {'0.1', '1', '10', '100'})
    legend(label);
    title(allSegments{s})
    
end



%% Plot Sensitivity for every stimulus contrast level

figure; hold all;
for s = [1,2]
    [val_tresh, idx_thresh] = min(abs(r(s).thisdPrime-threshold));
    sensitivity = 1./(expParams.contrastLevels(idx_thresh));
    
    
    plot([0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 6.5, 8],sensitivity);
    
end

set(gca,'XScale','log'); xlim([.1 100]); %ylim([1 1000]);
legend('Inner', 'Outer');
xlabel('Spatial frequency (cpd)')
ylabel('Contrast Sensitivity')
set(gca, 'XTick', [0.1, 1, 10, 100], 'XTickLabel', {'0.1', '1', '10', '100'})
set(gca, 'YTick', [1, 10, 100, 1000], 'YTickLabel', {'1', '10', '100', '1000'})


title('CSF eccen 10')
box off;

