% Banks et al. 1991

% Scene

% Horizontal Gabor patches
% Presented foveally or along the horizontal meridian

% Space average luminance - 762 cd/m2

% emmetrope and 0.75-D myope

% 1.5 mm artificial pupil

% Retinal illuminance was 1348 Td.

% 2 100 ms intervals of abrupt onset and offset

% 2AFC detection task, gabor or not.

% contrast varied

% ocular media transmittance, OTF - optics diffraction limited - identical
% at all eccentricities

% Curcio PR data


% SF = 1, then target size is 1.20 cycles per 2SD
% So 1 SD will be 0.6

%% Specify experiment parameters

nTrials         = 1;        % Number of trials per stimulus condition
contrast_levels = 1;        % Contrast levels
% eyemovements    = [1 1 0]';  % Which type of eye movments
eccentricities  = 10; %[0 2 5 10 20 40]; %4.5;
spatFreq          = 1; % [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26];
% polarangles     = 0;
% defocuslevels   = 0;         % units??  [0 0.5 1 1.5 2]
% verbose         = true;

% Temporal properties of one trial
tStep            = 0.01;                % Time step for optical image sequence (seconds)
tsamples         = (0:tStep:0.20);      % seconds
% timesd           = .01;              % sd of temporal Gaussian window

% Scene field of view
sparams.fov       = 0.5;   % scene field of view in degrees (diameter)
fov2deg           = 1/sparams.fov;

% Gabor parameters
tparams          = harmonicP;                   % Standard Gabor
tparams.ang       = pi/2;                % Gabor orientation (radians) - question: what is 0??
tparams.freq      = spatFreq*fov2deg;           % Spatial frequency (cycles/FOV)
tparams.contrast  = 1;                   % Presumably michelson, [0 1]
tparams.GaborFlag = 0.6*fov2deg;         % Gaussian window

% Make a default Optical Image Sequence. We need some of the
% parameters from this. It will be overwritten as we change the
% contrast and/or optics.

% The matched, zero contrast, harmonic parameters
tparams(1) = tparams;
tparams(2) = tparams(1); tparams(1).contrast = 0;

% And then we make a Gaussian temporal modulation that brings the stimulus
% on and off
tseries = zeros(length(tsamples),1);
tseries(5:15) = 1;
% stimWeights = ieScale(fspecial('gaussian',[1,200],50), 0, 1);

% add optics
oi = oiCreate('diffraction');

% The two harmonics are 'blended', which means at each moment in time we
% have a weighted sum of the two where the weights sum to 1.
[ois, scenes] = oisCreate('harmonic','blend',tseries, ...
    'testParameters',tparams,...
    'sceneParameters',sparams, ...
    'oi', oi, ...
    'meanluminance', 762);

% Visualize illuminance
ois.visualize('movie illuminance')

% Now, show the time series of weights
ois.visualize('weights');


%%

mosaicParams = struct(...
    'name', 'the hex mosaic', ...
    'resamplingFactor', 9, ...                      % Sets underlying pixel spacing; controls the rectangular sampling of the hex mosaic grid
    'eccBasedConeDensity', true, ...                % Whether to have an eccentricity based, spatially - varying density
    'sConeMinDistanceFactor', 3.0, ...              % Min distance between neighboring S-cones = f * local cone separation - used to make the S-cone lattice semi-regular
    'sConeFreeRadiusMicrons', 0.15*300, ...         % Radius of S-cone free retina, in microns (here set to 0.15 deg).
    'spatialDensity', [0 6/10 3/10 1/10]...         % With a LMS density of of 6:3:1
    );

quality.tolerance1 = 0.5;                           % larger than default tolerances to speed-up computation. For production work, either do not set, or set to equal or lower than 0.01
quality.tolerance2 = 0.05;                          % larger than default tolerances to speed-up computation, For production work, either do not set, or set to equal or lower than 0.001
quality.marginF = [];                               % How much larger lattice to generate so as to minimize artifacts in cone spacing near the edges. If empty, a dynamic adjustment of margin is done for mosaics < 1.0 degs

% Set import/export options
saveMosaic = false;                                 % whether to save the mosaic
loadMosaic = false;                                 % whether to load a previously saved mosaic
saveMosaicPDF = false;                              % whether to save a PDF of the mosaic

% Set FOVs examined
mosaicParams.fovDegs = sparams.fov;                 % mosaic FOV

theHexMosaic = coneMosaicHex(mosaicParams.resamplingFactor, ...
    'name', mosaicParams.name, ...
    'fovDegs', mosaicParams.fovDegs, ...
    'eccBasedConeDensity', mosaicParams.eccBasedConeDensity, ...
    'sConeMinDistanceFactor', mosaicParams.sConeMinDistanceFactor, ...
    'sConeFreeRadiusMicrons', mosaicParams.sConeFreeRadiusMicrons, ...
    'spatialDensity', mosaicParams.spatialDensity, ...
    'latticeAdjustmentPositionalToleranceF', quality.tolerance1, ...
    'latticeAdjustmentDelaunayToleranceF', quality.tolerance2, ...
    'marginF', quality.marginF ...
    );


% theHexMosaic.displayInfo();
% 
% % Visualize the mosaic, showing both the light collecting area (inner segment) and the geometric area
% visualizedAperture = 'lightCollectingArea'; % choose between 'both', 'lightCollectingArea', 'geometricArea'
% theHexMosaic.visualizeGrid(...
%     'visualizedConeAperture', visualizedAperture, ...
%     'apertureShape', 'disks', ...
%     'panelPosition', [1 1], 'generateNewFigure', true);


isomerizationsHex = theHexMosaic.compute(ois{1},'currentFlag',false);
theHexMosaic.window;





%%
spatialFrequencies  = [.25 .4 .65 1 1.6 2.6 4 6.5 8 10 16 26]; % cycles per degree

eccentricities      = [0 2 5 10 20 40];


[sf, eccen] = meshgrid(spatialFrequencies, eccentricities);
numCycles   = [...
    
NaN .8 1.05 1.2 1.35 1.5 1.75 1.9 NaN 2.03 2.05 2.05; ...   Fovea

NaN NaN 1 1.3 1.7 2.4 3.2 3.21 NaN 3.79 3.99 4.22; ...      2 deg

1 1.21 1.51 1.71 2.2 2.5 3.4 4.19 NaN 4.8 NaN NaN; ...      5 deg

0.8 1 1.49 2 2.39 3.02 3.99 4.5 4.5 NaN NaN NaN; ...        10 deg

1.7 1.7 1.81 1.90 2 2 2 NaN NaN NaN NaN NaN;  ...           20 deg

1.7 1.7 1.7 1.7 1.7 NaN NaN NaN NaN NaN NaN NaN; ...        40 deg

];



cycleLength = 1./sf;          % degrees per cycle

sd = cycleLength .* numCycles/2;



%%