% Visualize results

%% s_visualizeResults.m

% Define path where results are
resultsDir = fullfile(banksRootPath, 'results');

% General parameters
whichObserver = 'ideal';
allSegments = {'inner', 'outer'};

% Load experiment parameters
expParams = loadExpParamsBanks1991;

r_inner = cell(length(expParams.eccen(:,1)),1);
r_outer = r_inner;

for ec = 1:length(expParams.eccen(:,1))
    %
    
    inner = load(fullfile(resultsDir, sprintf('dPrime_%s_inner_eccen%d.mat', whichObserver, expParams.eccen(ec,1))));
    outer = load(fullfile(resultsDir, sprintf('dPrime_%s_outer_eccen%d.mat', whichObserver, expParams.eccen(ec,1))));
    r_inner{ec} =  inner.thisdPrime;
    r_outer{ec} =  outer.thisdPrime;
end






%%Plot contrast sensitivity functions for each eccentricity and cone segment
threshold = 1.36/sqrt(2);

for ii = 1:12, labelSF{ii} = sprintf('%1.3f',expParams.sf(1,ii)); end


cmap = jet(12);
allData = {r_inner,r_outer};

sensitivity = NaN(2,length(expParams.eccen(:,1)), length(expParams.sf(1,:)));

for s = [1,2]
    figure(s); set(gcf, 'Color', 'w'); clf;
    dprime = allData{s};
    
    for ec = 1:length(expParams.eccen(:,1))
        
        idx = find(~isnan(expParams.sd(ec,:)'));
        
        
        subplot(3,2,ec); hold all;
        for ii = idx'
            plot(expParams.contrastLevels*100, dprime{ec}(:,ii), 'Color', cmap(ii,:,:), 'LineWidth',2)
            [val_tresh, idx_thresh] = min(abs(dprime{ec}(:,ii) -threshold));
            sensitivity(s, ec, ii) = 1./(expParams.contrastLevels(idx_thresh));
        end
        
        
        
        plot([0.1 100], [threshold, threshold],'k--')
        % plot([0.1 100], [1.36, 1.36],'r--')
        xlim([0.1 100]); ylim([0.1 100]);
        xlabel('Contrast (Michelson)')
        ylabel('Contrast Sensitivity (d prime)'); title('Figure 3')
        set(gca, 'XScale', 'log', 'YScale', 'log')
        set(gca, 'XTick', [0.1, 1, 10, 100], 'XTickLabel', {'0.1', '1', '10', '100'})
        set(gca, 'YTick', [0.1, 1, 10, 100], 'YTickLabel', {'0.1', '1', '10', '100'})
        legend({labelSF{idx}, 'threshold'}, 'Location', 'Best');
        title(sprintf('%s: eccen %d', allSegments{s}, expParams.eccen(ec,1)))
    end
    
end


%% Plot Sensitivity for every stimulus contrast level

for ii = 1:6, labelEccen{ii} = sprintf('%1.3f',expParams.eccen(ii,1)); end


for s = [1,2]
    figure; hold all;
    
    for ec = 1:length(expParams.eccen(:,1))
        
        plot(expParams.sf(1,:),squeeze(sensitivity(s,ec,:)), 'LineWidth',2);        
    end
    
    set(gca,'XScale','log'); xlim([.1 100]); ylim([.1 1000]);
    xlabel('Spatial frequency (cpd)')
    ylabel('Contrast Sensitivity')
    set(gca, 'XTick', [0.1, 1, 10, 100], 'XTickLabel', {'0.1', '1', '10', '100'}, 'TickDir', 'out')
    set(gca, 'YTick', [0.1, 1, 10, 100, 1000], 'YTickLabel', {'0.1','1', '10', '100', '1000'}, 'TickDir', 'out')
    
    title(allSegments{s})
    box off;
    legend(labelEccen)
    
end




