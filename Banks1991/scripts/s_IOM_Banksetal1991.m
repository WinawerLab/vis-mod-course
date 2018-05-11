%% s_IOM_Banksetal1991

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
%   - Check inter cone spacing (Fig 1) - Check Retinal coverage   (Fig 2) -
%   Check Proportion absorbed (Fig 2)


% Get size of stimuli in degrees per 1 SD from table 1
[sd, sf, eccenStim, cycleLength] = table1SizetoSD;

verbose = false;
%% 0. Specify experiment parameters

nTrials         = 1;        % Number of trials per stimulus condition
contrastLevels  = [0.01:0.01:0.1, 0.2:0.1:1];        % Contrast levels
eccentricities  = 10;       % [0 2 5 10 20 40];
spatFreq        = 0.25;        % [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26];
% polarangles     = 0; % Horizontal meridian defocuslevels   = 0; % units??
% eyemovements    = [0 0 0]';  % No eye movments verbose         = true;


% Temporal properties of one trial
tStep            = 0.01;                % Time step for optical image sequence (seconds)
tSamples         = (0:tStep:0.20);      % seconds
% timesd           = .01;               % sd of temporal Gaussian window

% Scene field of view
sparams.fov       = 2;   % scene field of view in degrees (diameter)
fov2deg           = 1/sparams.fov;
sparams.distance  = 0.57;  % meters

% Gabor parameters
tparams           = harmonicP;           % Function to get standard Gabor
tparams.ang       = pi/2;                % Gabor orientation (radians) - for now, horizontal
tparams.freq      = spatFreq*fov2deg;    % Spatial frequency (cycles/FOV)
tparams.GaborFlag = 1.6/2;               % Gaussian window (1sd?)

% Define gabors as test params, and add another one for the blank stimulus
tparams(1) = tparams;
tparams(1).contrast  = 0;                   % Presumably michelson, [0 1]
tparams(2) = tparams(1);

% Make a Gaussian temporal modulation that brings the stimulus on and off
tseries = zeros(length(tSamples),1);
tseries(5:15) = 1; % HACK
% stimWeights = ieScale(fspecial('gaussian',[1,200],50), 0, 1);

% Set Mosaic import/export options
saveMosaic = false;                                 % whether to save the mosaic
loadMosaic = false;                                 % whether to load a previously saved mosaic
saveMosaicPDF = false;                              % whether to save a PDF of the mosaic

quality.tolerance1 = 0.5;                           % larger than default tolerances to speed-up computation. For production work, either do not set, or set to equal or lower than 0.01
quality.tolerance2 = 0.05;                          % larger than default tolerances to speed-up computation, For production work, either do not set, or set to equal or lower than 0.001
quality.marginF = [];                               % How much larger lattice to generate so as to minimize artifacts in cone spacing near the edges. If empty, a dynamic adjustment of margin is done for mosaics < 1.0 degs


%% 1. Create optics
% Add optics (Question: how to add pupil size?)
wvf = wvfCreate('calc wave', [400:10:700]);
wvf = wvfSet(wvf,'measured pupil size', 1.5);
wvf = wvfSet(wvf,'calc pupil size', 1.5);
wvf = wvfComputePSF(wvf);

% plot PSF
% wvfPlot(wvf, '1d psf angle', 'deg', 550)
% wvfPlot(wvf, 'otf', 'deg', 550)

oi = wvf2oi(wvf);
oi = opticsSet(oi, 'model', 'diffraction limited');

% Plot ocular transmittance function
% figure; plot(oi.optics.lens.wave, oi.optics.lens.transmittance)
title('Ocular transmittance function')
% % Request pupil diameter
% p = opticsGet(oi.optics,'pupil diameter','mm');
% 
% if verbose
%     % Plot OTF
%     oiPlot(oi, 'OTF', 'wavelength', 550)
%     
%     % Plot Point spread function
%     oiPlot(oi, 'PSF', 'wavelength', 550)
%     
%     % Plot line spread function for multiple wave lengths
%     oiPlot(oi,'ls wavelength');
%     title(sprintf('F/# = %.0d',opticsGet(oi.optics,'f number')))
% end

%% 2. Create OI Sequence

% The two harmonics are 'blended', which means at each moment in time we
% have a weighted sum of the two where the weights sum to 1.
[ois(1), scene1] = oisCreate('harmonic','blend',tseries, ...
    'testParameters',tparams,...
    'sceneParameters',sparams, ...
    'oi', oi, ...
    'meanluminance', 762);

ois(1).visualize('movieilluminance');

ieAddObject(scene1{1});
ieAddObject(scene1{2});
% sceneWindow;

%% 3. Create Cone Mosaic

whichEye = 'left';
deg2m    = 0.3 * 0.001; % 3 deg per mm, .001 mm per meter
polAng   = 0;




%%
alpha_absorptions = cell(length(eccentricities));
beta_absorptions = alpha_absorptions;

for eccen = eccentricities
    
    % Compute x,y position in m of center of retinal patch from ecc and angle
    [x, y] = pol2cart(polAng, eccen);
    x = x * deg2m;  y = y * deg2m;
    
    regMosaicParams = struct( ...
        'eccentricity', eccen, ...
        'polarAngle', polAng, ... Right horizontal meridian
        'cmFOV', sparams.fov);
    
    cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
    
    % Set the field of view (degrees)
    cMosaic.setSizeToFOV(regMosaicParams.cmFOV);
    
    % Add photon noise
    cMosaic.noiseFlag = 'frozen'; % 'random' 'frozen' 'none'
    
    % Not sure why these have to match, but there is a bug if they don't.
    cMosaic.integrationTime = ois(1).timeStep;
    
    % There are no eyemovements, but I think you need to have emPaths defined in
    % order to get time varying absorption rates (because it's an oisSequence)
    regMosaicParams.em        = emCreate;
    regMosaicParams.em.emFlag =  [0 0 0]';
    emPaths  = cMosaic.emGenSequence(length(tSamples), 'nTrials', nTrials, ...
        'em', regMosaicParams.em);
    cMosaic.emPositions = emPaths;
   
    % implement th inner segment aperture to correct for proportion covered
    propCovered = getBanks1991ConeCoverage(eccen);

    cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered;
    cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered;

    
    for c = contrastLevels
        
        % recompute stim for particular contrast
        tparams(2).contrast = c;
        
        [ois(2), scene2] = oisCreate('harmonic','blend',tseries, ...
            'testParameters',tparams,...
            'sceneParameters',sparams, ...
            'oi', oi, ...
            'meanluminance', 762);
     
        if verbose
            % Visualize illuminance
            ois(1).visualize('movie illuminance')
            ois(2).visualize('movie illuminance')
            
            % Now, show the time series of weighting the Gabor and blank stimulus
            ois(2).visualize('weights');
        end
        
        % Compute absorptions
        alpha_absorptions{eccen==eccentricities}(:,:,:,:,c==contrastLevels) = cMosaic.compute(ois(1), 'currentFlag', false, 'emPaths', emPaths);
        beta_absorptions{eccen==eccentricities}(:,:,:,:,c==contrastLevels) = cMosaic.compute(ois(2), 'currentFlag', false, 'emPaths', emPaths);
        
    end
end


%% Visualize cone mosaic absorptions
% ieAddObject(scene2{1})
% ieAddObject(scene2{2})
% sceneWindow;
% cMosaic.window;



%% Calculate d prime

% Likelihood
dPrimeFunction1 = @(alpha, beta)  (sum( (beta(:)-alpha(:)) .* log(beta(:)./alpha(:))  ) ./ ...
                                sqrt(0.5* sum( (beta(:)+alpha(:)) .* log( (beta(:)./alpha(:))).^2 )) ); 

dPrimeFunction2 = @(alpha, beta) ((nanmean(beta(:)-alpha(:)))./sqrt(nanmean(beta(:))));

% Intensity
dPrimeFunction3 = @(alpha, beta) (1.36*sqrt(nanmean(beta(:))));

N = @(alpha, beta) nansum((beta(:)+alpha(:))/2);
deltaN = @(alpha, beta) nansum(beta(:)-alpha(:));

dPrimeFunction4 = @(alpha, beta) deltaN(alpha, beta)/N(alpha, beta);

d_prime1 = zeros(length(eccentricities), length(contrastLevels));
% d_prime2 = d_prime1; d_prime3 = d_prime1; d_prime4 = d_prime1;

for eccen = eccentricities
    for c=contrastLevels
        
        this_alpha = squeeze(mean(alpha_absorptions{eccen==eccentricities}(1,:,:,:,c==contrastLevels),4));
        this_beta = squeeze(mean(beta_absorptions{eccen==eccentricities}(1,:,:,:,c==contrastLevels),4));
        
        numerator = nansum( (this_beta(:) - this_alpha(:)) .* log(this_beta(:)./this_alpha(:)) );
        denumerator = sqrt(0.5* nansum( (this_beta(:) + this_alpha(:)) .* log(this_beta(:)./this_alpha(:)).^2 ));
        check_d(eccen==eccentricities,c==contrastLevels) = numerator/denumerator;
        
        d_prime1(eccen==eccentricities,c==contrastLevels) = dPrimeFunction1(this_alpha,this_beta);
%         d_prime2(eccen==eccentricities,c==contrastLevels) = dPrimeFunction2(this_alpha,this_beta);
%         d_prime3(eccen==eccentricities,c==contrastLevels) = dPrimeFunction3(this_alpha,this_beta);
%         d_prime4(eccen==eccentricities,c==contrastLevels) = dPrimeFunction4(this_alpha,this_beta);
    end
end

threshold = 1.36*sqrt(2);

figure(3); clf;
plot(contrastLevels, d_prime1); 
% plot(log10(contrastLevels), log10(d_prime1)); 
hold on; plot([0 1], ([threshold, threshold]),'k--')
% hold on; plot([-2 1], log10([1.36, 1.36]),'k--')
xlabel('Contrast (Michelson)')
ylabel('D prime'); title('Calculation using Log likelihood')
% 
% subplot(222); plot(contrastLevels,check_d);
% hold on; plot([0 1], [threshold, threshold],'k--')
% xlabel('Contrast (Michelson)')
% ylabel('D prime'); title('Calculation check')
% 
% subplot(223); plot(contrastLevels, d_prime3);
% xlabel('Contrast (Michelson)')
% ylabel('D prime'); title('Calculation using 1.36*sqrt(N)')
% 
% subplot(224); plot(contrastLevels, d_prime4);
% xlabel('Contrast (Michelson)')
% ylabel('D prime'); title('Calculation using delta N / N')



return

%% Hex Mosaic
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

