%% s_visualizeResults.m

% Script to visualize Figure 3 and 5 (contrast sensitivity functions for 
% different spatial frequencies and eccentricities for ideal observer)
% from the published article by Banks et al. (1991).

% Full citation:
% Peripheral spatial vision: limits imposed by optics, photoreceptors, and
% receptor pooling. (1991) Banks, M.S., Sekuler, S.B., & Anderson, S.J.
% Journal of the Optical Society of America A (JOSAA). Vol 8, No 11,
% 1775-1787.

% DESCRIPTION:
% This script recompute Figure 3 and 5 as follows:
%   0. Set up general parameters and experimental parameters:
%       - Path where results live (computed by s_mainAnalysis.m)
%       - Load stimulus parameters, like size, spatial frequency,
%       eccentricity
%   1. Plot d-prime against stimulus contrast, for every spatial frequency,
%      for every eccentricity (0, 2, 5, 10, 20 degrees).
%   2. Get contrast threshold for each spatial frequency, for each
%      eccentricity. Plot the threshold's against spatial frequency.


% To get results use s_mainAnalysis.m (This script requires ISETBio)

% DEPENDENCIES:
%   Results computed by s_mainAnalysis.m

% AUTHOR:
%   Eline Kupers, 2018, NYU



%% 0. GENERAL AND EXPERIMENT PARAMETERS 

% Define path where results are
resultsDir = fullfile(banksRootPath, 'results');
figureDir = fullfile(banksRootPath, 'figures');
saveFigures = true;

% General parameters
whichObserver = 'ideal';
allSegments = {'inner','outer'};

% Load experiment parameters
expParams = loadExpParamsBanks1991;

% Preallocate space for dPrime
r_inner = cell(length(expParams.eccen(:,1)),1);
r_outer = r_inner;

% For every eccentricity..
for ec = 1:length(expParams.eccen(:,1))
    
    % Load and store inner and outer segment results
    inner = load(fullfile(resultsDir, sprintf('dPrime_%s_inner_eccen%d.mat', whichObserver, expParams.eccen(ec,1))));
    outer = load(fullfile(resultsDir, sprintf('dPrime_%s_outer_eccen%d.mat', whichObserver, expParams.eccen(ec,1))));
    r_inner{ec} =  inner.thisdPrime;
    r_outer{ec} =  outer.thisdPrime;
    
end



%% 1. PLOT CONTRAST LEVELS as a function of ECCENTRICITY and CONE SEGMENT (Fig 3)

% Ideal observer model threshold (Based on Geisler (1987))
threshold = 1.36/sqrt(2);

% Get labels 
for ii = 1:12, labelSF{ii} = sprintf('%1.2f',expParams.sf(1,ii)); end

% Get color map
cmap = jet(size(expParams.sf,2));

% Get all data
allData = {r_inner,r_outer};

% Preallocate space for contrast sensitivity
sensitivity = NaN(2,length(expParams.eccen(:,1)), length(expParams.sf(1,:)));

% For each cone segment..
for s = [1,2]
    
    % Create a figure
    figure(s); set(gcf, 'Color', 'w', 'Position', [28, 148, 1072, 1193]); clf;
    
    % Select d-prime data
    dprime = allData{s};
    
    % For every eccentricity
    for ec = 1:length(expParams.eccen(:,1))
        
        % Check for NaNs
        idx = find(~isnan(expParams.sd(ec,:)'));
        
        % Plot stimulus contrast against d-prime
        subplot(3,2,ec); hold all;
        for ii = idx'
            plot(expParams.contrastLevels*100, dprime{ec}(:,ii), 'Color', cmap(ii,:,:), 'LineWidth',2)
            [val_tresh, idx_thresh] = min(abs(dprime{ec}(:,ii) -threshold));
            sensitivity(s, ec, ii) = 1./(expParams.contrastLevels(idx_thresh));
        end
        
        
        % Add labels & make plot pretty
        plot([0.1 100], [threshold, threshold],'k--')
        % plot([0.1 100], [1.36, 1.36],'r--')
        xlim([0.1 100]); ylim([0.1 100]);
        xlabel('Stimulus contrast (%)')
        ylabel('Sensitivity (d-prime)'); title('Figure 3')
        set(gca, 'XScale', 'log', 'YScale', 'log')
        set(gca, 'XTick', [0.1, 1, 10, 100], 'XTickLabel', {'0.1', '1', '10', '100'},'TickDir', 'out', 'FontSize', 20)
        set(gca, 'YTick', [0.1, 1, 10, 100], 'YTickLabel', {'0.1', '1', '10', '100'},'TickDir', 'out', 'FontSize', 20)
        legend({labelSF{idx}, 'threshold'}, 'FontSize', 10,'Location', 'Best'); legend boxoff
        title(sprintf('%s: eccen %d', allSegments{s}, expParams.eccen(ec,1)))
    end
    
    if saveFigures
        print(fullfile(figureDir, sprintf('Fig3_dPrime_%s_%s', whichObserver,  allSegments{s})),'-dpng');
    end
end



%% 2. PLOT CONTRAST SENSITIVITY as a function of SPATIAL FREQUENCY (Fig 5)

% Get labels
for ii = 1:6, labelEccen{ii} = sprintf('%d',expParams.eccen(ii,1)); end

% Get markers
M = {'^', 'o', 's', 's', '^', 'o'};

edgeColors = jet(length(expParams.eccen(:,1)));
faceColors = edgeColors;
faceColors(1,:) = [1 1 1];
faceColors(4,:) = [1 1 1];
faceColors(6,:) = [1 1 1];

% Set up figure
figure(3); clf; set(gcf, 'Color', 'w', 'Position', [62, 611, 1217, 734]);

% For each cone segment
for s = [1,2]
    
    % Get a subplot
    subplot(1,2,s); hold all;
    
    % For each eccentricity
    for ec = 1:length(expParams.eccen(:,1))        
        plot(expParams.sf(1,:),squeeze(sensitivity(s,ec,:)), ...
            'Color', [edgeColors(ec,:)], ...
            'Marker', M{ec}, 'MarkerSize', 20, ...
            'MarkerEdgeColor', edgeColors(ec,:), ... 
            'MarkerFaceColor', faceColors(ec,:), ...
            'LineWidth',2);        
    end
    
    % Add labels & make plot pretty
    set(gca,'XScale','log'); xlim([0.1 100]); ylim([1 1000]);
    xlabel('Spatial frequency (cpd)')
    ylabel('Contrast Sensitivity (%)')
    set(gca, 'XTick', [0.1, 1, 10, 100], 'XTickLabel', {'0.1', '1', '10', '100'}, 'TickDir', 'out', 'FontSize', 20)
    set(gca, 'YTick', [1, 100, 1000], 'YTickLabel', {'1', '100', '1000'}, 'TickDir', 'out', 'FontSize', 20)    
    title([allSegments{s} ' segment aperture'])
    box off;
    legend(labelEccen, 'FontSize',20, 'Location', 'Best'); legend boxoff
    
end

if saveFigures
    print(fullfile(figureDir, sprintf('Fig5_SpatFreqLimits_%s', whichObserver)),'-dpng');
end


