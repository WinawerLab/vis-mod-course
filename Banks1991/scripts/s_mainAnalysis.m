%% s_mainAnalysis.m

% Script to create an ideal observer model based on the article by Banks et
% al. 1991

% This toolbox requires ISETBio

% By Eline Kupers, 2018

%% Description & details


% SCENE DETAILS:

% - The stimuli are horizontal Gabor patches, presented foveally or along
% the horizontal meridian 2AFC detection task: One Gabor presented abruptly
% for 100 ms intervals and one blank screen. Task was to say which interval
% there was a Gabor - Average luminance = 762 cd/m2 - Retinal illuminance
% was 1348 Td. - Contrast levels were varied between 1-100%

% OPTICS DETAILS
% - Subjects are one emmetrope and one 0.75-D myope - 1.5 mm
% artificial pupil
% - OTF was optics diffraction limited and identical at
% all eccentricities

% NOTE: Banks et al only takes Pupil size, ocular media transmittance and
% photoreceptor properties of the human eye at 0, 2, 5, 10, 20, 40 deg
% eccentricity (horizontal VF, nasal retina).

% PHOTORECEPTOR DETAILS - Curcio PR data was used to model inner cone
% segments

% NEXT STEPS / TO DO:
% - Change pupil size to 1.5 mm
% - Figure out trade off mosaic resampling versus computing time
% - Related: how do we match scene FoV and stim size with cone mosaic FoV?
% - How to calculate d'? Do we sum or average over time?


% Recompute Figure 3, as follows: For a given eccentricity (say 10 degrees)
% and spatial frequency (say 8 cpd)
%  - test target sizes 1-10 (use convertion of table 1 to get Gaussian SD)
%  - then vary contrast levels
%  - compute cone absorptions for blank and Gabor stimuli
%  - Use a linear classifier on the absorption levels
%  - Use the classifier performance to plot the psychometric functions
%  - Take the contrast sensitivity threshold (in d prime?) at 75%, plot
%  against target size for each SF and Eccen.


% For a given eccentricity:
%   - Check inter cone spacing (Fig 1) - We currently use the inner segment
%   - Proportion absorbed (Fig 2) -- See propCovered = getBanks1991ConeCoverage(thisEccen);



%% 0. General parameters
verbose = false;
deg2m    = 0.3 * 0.001; % 3 deg per mm, .001 mm per meter

whichObserver = 'ideal'; % choose from 'ideal' or 'human'
segment       = 'inner'; % choose from 'inner' or 'outer' cone segment
saveResults   = true;
%% 1. Specify experiment parameters

% Load experiment parameters
expParams = loadExpParamsBanks1991;

% Set a dummy contrast level to create stimulus test params
thisContrast = expParams.contrastLevels(1);

for thisEccen = 10;%expParams.eccen(:,1)'
    % thisEccen    = 10; % Choose from 0, 2, 5, 10, 20, 40
     
    for thisSpatFreq = [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 6.5, 8]; %expParams.sf(1,:)
        % thisSpatFreq = 0.25;  % Choose from [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 6.5, 8, 10, 16, 26]
        
        % Find the indices of corresponding target size
        idx = [find(thisEccen == expParams.eccen(:,1)), find(thisSpatFreq == expParams.sf(1,:))];
        
        % Get the Gaussian window of the stimulus target
        thisTargetSize = expParams.sd(idx(1),idx(2));
        
        if ~isnan(thisTargetSize)
            % Load stimulus params
            [tparams, sparams, tseries] = getStimulusParams(thisContrast, thisTargetSize, thisSpatFreq);
        
        
            %% 2. Create optics
            oi = getOptics(whichObserver, verbose);

            %% 3. Create OI Sequence

            % The two harmonics are 'blended', which means at each moment in time we
            % have a weighted sum of the two where the weights sum to 1.
            [ois(1), scene1] = oisCreate('harmonic','blend',tseries, ...
                'testParameters',tparams{1},...
                'sceneParameters',sparams, ...
                'oi', oi, ...
                'meanluminance', 762);

            %% 4. Create Cone Mosaic

            % Preallocate cells for absorptions
%             alphaAbsorptions = zeros(length(expParams.contrastLevels));
%             betaAbsorptions = []; %alpha_absorptions;


            cMosaic = getConeMosaic(expParams, thisEccen, deg2m, sparams, ois, whichObserver, segment, verbose);

            for c = expParams.contrastLevels

                % recompute stim for particular contrast
                tparams{2}(2).contrast = c;

                [ois(2), scene2] = oisCreate('harmonic','blend',tseries, ...
                    'testParameters',tparams{2},...
                    'sceneParameters',sparams, ...
                    'oi', oi, ...
                    'meanluminance', 762);

                if verbose
           
                    % visualize the OIS
                    ois(2).visualize('movieilluminance');

                    % visualize the scene
                    ieAddObject(scene2{1});
                    ieAddObject(scene2{2});
                    sceneWindow;

                    % Now, show the time series of weighting the Gabor and blank stimulus
%                     ois(2).visualize('weights');
                end

                % Compute absorptions 6D
                % (contrast x SF x trials x cols x rows x time points)
                alphaAbsorptions(c==expParams.contrastLevels, thisSpatFreq==expParams.sf(1,:),:,:,:,:) = cMosaic.compute(ois(1), 'currentFlag', false);
                betaAbsorptions(c==expParams.contrastLevels, thisSpatFreq==expParams.sf(1,:),:,:,:,:) = cMosaic.compute(ois(2), 'currentFlag', false);

            end
            
        end
    end
end

% Visualize cone mosaic absorptions
if verbose; cMosaic.window; end

%% 5. Calculate sensitivity (d-prime)

% Likelihood
dPrime = dPrimeFunction;

thisdPrime = [];
for c = 1:length(expParams.contrastLevels)
    for sf=1:9

        this_alpha = squeeze(mean(alphaAbsorptions(c,sf, 1,:,:,:),6));
        this_beta = squeeze(mean(betaAbsorptions(c,sf, 1,:,:,:),6));

        thisdPrime(c, sf) = dPrime(this_alpha,this_beta);

    end
end


if saveResults
    saveDir = fullfile(banksRootPath, 'results'); 
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    save(fullfile(saveDir, sprintf('dPrime_%s_%s_eccen%d.mat', whichObserver, segment, thisEccen)), 'thisdPrime');   
end
    




return

%% Hex Mosaic

% Set Hexagonal Mosaic import/export options
saveMosaic = false;                                 % whether to save the mosaic
loadMosaic = false;                                 % whether to load a previously saved mosaic
saveMosaicPDF = false;                              % whether to save a PDF of the mosaic

quality.tolerance1 = 0.5;                           % larger than default tolerances to speed-up computation. For production work, either do not set, or set to equal or lower than 0.01
quality.tolerance2 = 0.05;                          % larger than default tolerances to speed-up computation, For production work, either do not set, or set to equal or lower than 0.001
quality.marginF = [];                               % How much larger lattice to generate so as to minimize artifacts in cone spacing near the edges. If empty, a dynamic adjustment of margin is done for mosaics < 1.0 degs



HexMosaicParams = struct(...
    'name', 'the hex mosaic', ...
    'resamplingFactor', 9, ...                      % Sets underlying pixel spacing; controls the rectangular sampling of the hex mosaic grid
    'eccBasedConeDensity', true, ...                % Whether to have an eccentricity based, spatially - varying density
    'sConeMinDistanceFactor', 3.0, ...              % Min distance between neighboring S-cones = f * local cone separation - used to make the S-cone lattice semi-regular
    'sConeFreeRadiusMicrons', 0.15*300, ...         % Radius of S-cone free retina, in microns (here set to 0.15 deg).
    'spatialDensity', [0 6/10 3/10 1/10]...         % With a LMS density of of 6:3:1
    );

% Set FOVs examined
HexMosaicParams.fovDegs = sparams.fov;                 % mosaic FOV

% Create the hexagonal mosaic
theHexMosaic = coneMosaicHex(HexMosaicParams.resamplingFactor, ...
    'name', HexMosaicParams.name, ...
    'fovDegs', HexMosaicParams.fovDegs, ...
    'eccBasedConeDensity', HexMosaicParams.eccBasedConeDensity, ...
    'sConeMinDistanceFactor', HexMosaicParams.sConeMinDistanceFactor, ...
    'sConeFreeRadiusMicrons', HexMosaicParams.sConeFreeRadiusMicrons, ...
    'spatialDensity', HexMosaicParams.spatialDensity, ...
    'latticeAdjustmentPositionalToleranceF', quality.tolerance1, ...
    'latticeAdjustmentDelaunayToleranceF', quality.tolerance2, ...
    'marginF', quality.marginF ...
    );


theHexMosaic.displayInfo();

% Visualize the mosaic, showing both the light collecting area (inner
% segment) and the geometric area
visualizedAperture = 'lightCollectingArea'; % choose between 'both', 'lightCollectingArea', 'geometricArea'
theHexMosaic.visualizeGrid(...
    'visualizedConeAperture', visualizedAperture, ...
    'apertureShape', 'disks', ...
    'panelPosition', [1 1], 'generateNewFigure', true);

%% 4. Compute the cone isomerizations / absorptions
isomerizationsHex = theHexMosaic.compute(ois,'currentFlag',false);


% Visualize
theHexMosaic.window;

