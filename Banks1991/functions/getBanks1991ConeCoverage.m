function propCovered = getBanks1991ConeCoverage(eccentricity, varargin)

% Syntax
%   propCovered = getBanks1991ConeCoverage(eccentricity, ['eccentricityUnits'], ['deg'])

% Get proportion covered by cones for different eccentricities (ratio), by fitting
% a straight line in semi-log space to data in Figure 2 of Banks et al.
% (1991): Peripheral spatial vision: limited imposed by optics,
% photoreceptors and receptor pooling. J. Opt. Soc Am.

% Inputs
% eccentricity:         integer defining the eccentricity of interest, can 
%                       only be integer: 0, 2.5, 5, 10, 20, or 40 deg
% eccentricityUnits:    [=optional], string defining unit of input
%                       eccentricity: 'deg' [default], 'm','mm','um') will 
%                       be converted to 'deg'.
% coneSegment:          [=optional], string defining whether to use inner
%                       or outer cone segment. 'inner' [default], or
%                       'outer'


% Outputs
% propCovered:          integer of proportion covered by the cones for the
%                       requested eccentricity
%
% Example:
%  propCovered = getBanks1991ConeCoverage(4)

%% Handle input arguments
p = inputParser;
p.KeepUnmatched = true;

% Required and added
p.addRequired('eccentricity',@isvector);
p.addParameter('eccentricityUnits','deg',@ischar);
p.addParameter('coneSegment','inner',@ischar);
p.addParameter('verbose',false);


% Parse
p.parse(eccentricity, varargin{:});

% Get corrent eccenctricity units
switch (p.Results.eccentricityUnits)
    
    case 'm'
        eccDeg = p.Results.eccentricity*(1000/0.3);
    case 'mm'
        eccDeg = p.Results.eccentricity*(1/0.3);
    case 'um'
        eccDeg = p.Results.eccentricity*(0.001/0.3);
    case 'deg'
        eccDeg = p.Results.eccentricity;
    otherwise
        error('(getBanks1991ConeSpacing): Unknonwn units for eccentricity specified');
end

% set manually:
allEccentricities = 0:0.5:50;
banksEccen = [0, 2, 5, 10, 20, 40];

switch p.Results.coneSegment
    
    % Data is inferred from figure, not actual from a publicly shared dataset    
    case 'inner'   
        banksData  = [1, 0.8, 0.5, 0.4,0.3, 0.25];
        
        % coverage = (exp(-1/40*allEccentricities));
        %   or with fitting function:
        f = fit(banksEccen', banksData','exp2');

        coverage = f.a*exp(f.b*allEccentricities) + f.c*exp(f.d*allEccentricities);
        coverage(1) = 1; %reset fovea to 1.

    case 'outer'
        banksData  = [0.25, 0.07, 0.045, 0.04, 0.025, 0.02];  
        
        f = fit(banksEccen', banksData','exp2');
        coverage = f.a*exp(f.b*allEccentricities) + f.c*exp(f.d*allEccentricities);

        
end

    
% Plot for debugging
if p.Results.verbose
    figure(101); clf; set(gcf, 'Color','w')
    plot(allEccentricities, coverage, 'o-'); set(gca, 'YScale', 'log', ...
        'XScale', 'linear', 'XLim', [0 50], 'YLim', [0.01 1], 'TickDir', 'out');
    hold all; scatter(banksEccen, banksData, 80','k');
    xlabel('Eccentricity (deg)'); ylabel('Normalized proportion absorpted');
    title(sprintf('Proportion absorpted as a function of eccentricity, using %s segment', p.Results.coneSegment))
    box off
end

idx = (eccDeg == allEccentricities);
if any(idx)
     propCovered = sqrt(coverage(idx));
else
    error('(getBanks1991ConeCoverage): No data is available for requested eccentricity');
end
    

return 