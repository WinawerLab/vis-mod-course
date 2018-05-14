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
%   - Check inter cone spacing (Fig 1)
%   - Check Retinal coverage   (Fig 2)
%   - Check Proportion absorbed (Fig 2)



%% 0. General parameters
verbose = false;
deg2m    = 0.3 * 0.001; % 3 deg per mm, .001 mm per meter

whichObserver = 'ideal'; % choose between 'ideal' or 'human'

dPrime = [];
dPrime2 = [];
dPrime4 = [];
Nalone = [];
%% 1. Specify experiment parameters

% Load experiment parameters
expParams = loadExpParamsBanks1991;

% What stim params to use
thisContrast = expParams.contrastLevels(1);

for thisEccen = expParams.eccen(:,1)'
    % thisEccen    = 10; % Choose from 0, 2, 5, 10, 20, 40
     
    for thisSpatFreq = expParams.sf(1,3:5)
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

            if verbose
                % visualize the OIS
                ois(1).visualize('movieilluminance');

                % visualize the scene
                ieAddObject(scene1{1});
                ieAddObject(scene1{2});
                sceneWindow;
            end

            %% 4. Create Cone Mosaic

            % Preallocate cells for absorptions
            alpha_absorptions = []; %zeros(length(expParams.contrastLevels));
            beta_absorptions = []; %alpha_absorptions;


            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(expParams.polarangle, thisEccen);
            x = x * deg2m;  y = y * deg2m;

            regMosaicParams = struct( ...
                'eccentricity', thisEccen, ...
                'polarAngle', expParams.polarangle, ... Right horizontal meridian
                'cmFOV', sparams.fov);

            cMosaic = coneMosaic('center', [x, y], 'whichEye', expParams.whichEye);

            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(regMosaicParams.cmFOV);

            % Add photon noise
            cMosaic.noiseFlag = 'random';


            % Not sure why these have to match, but there is a bug if they don't.
            cMosaic.integrationTime = ois(1).timeStep;

            % There are no eyemovements, but I think you need to have emPaths defined in
            % order to get time varying absorption rates (because it's an oisSequence)
            regMosaicParams.em        = emCreate;
            regMosaicParams.em.emFlag =  [0 0 0]';
            emPaths  = cMosaic.emGenSequence(ois(1).length, 'nTrials', expParams.nTrials, ...
                'em', regMosaicParams.em);
            cMosaic.emPositions = emPaths;

            % implement th inner segment aperture to correct for proportion covered
            propCovered = getBanks1991ConeCoverage(thisEccen);

            cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered;
            cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered;


            for c = expParams.contrastLevels(1)

                % recompute stim for particular contrast
                tparams{2}(2).contrast = c;

                [ois(2), scene2] = oisCreate('harmonic','blend',tseries, ...
                    'testParameters',tparams{2},...
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
                alpha_absorptions(:,:,:,:,c==expParams.contrastLevels,thisSpatFreq==expParams.sf(1,:), thisEccen==expParams.eccen(:,1)') = cMosaic.compute(ois(1), 'currentFlag', false, 'emPaths', emPaths);
                beta_absorptions(:,:,:,:,c==expParams.contrastLevels,thisSpatFreq==expParams.sf(1,:), thisEccen==expParams.eccen(:,1)') = cMosaic.compute(ois(2), 'currentFlag', false, 'emPaths', emPaths);

            end
            
        end
    end
end   
    % Visualize cone mosaic absorptions
    if verbose; cMosaic.window; end
    
    %% 5. Calculate sensitivity (d-prime)
    
    % Likelihood
    dPrimeFunction1 = @(alpha, beta)  (sum( (beta(:)-alpha(:)) .* log(beta(:)./alpha(:))  ) ./ ...
        sqrt(0.5* sum( (beta(:)+alpha(:)) .* log( (beta(:)./alpha(:))).^2 )) );
    
    d_prime1 = zeros(length(expParams.contrastLevels),length(expParams.eccen(:,1)'), length(expParams.sf(1,:)));
    
    for ec = expParams.eccen(:,1)'
        for sf=expParams.sf(1,:)
            
            this_alpha = squeeze(mean(alpha_absorptions(1,:,:,:, 1,ec==expParams.eccen(:,1)', sf==expParams.sf(1,:)),4));
            this_beta = squeeze(mean(beta_absorptions(1,:,:,:,1, ec==expParams.eccen(:,1)', sf==expParams.sf(1,:)),4));
            
            d_prime1(c==expParams.contrastLevels, ec==expParams.eccen(:,1)', sf==expParams.sf(1,:)) = dPrimeFunction1(this_alpha,this_beta);
        end
    end
    
    dPrime = [dPrime, d_prime1];
    


% dPrime = dPrime(1:28);
% threshold = 1.36*sqrt(2);


%%

% figure; set(gcf, 'Color', 'w'); clf; hold on;
% cmap = jet(9);
% for ii = 1:9
%     plot(expParams.contrastLevels*100, dPrime(:,ii), 'Color', cmap(ii,:,:))
% end%, 'LineWidth',4);
% % plot(log10(contrastLevels), log10(d_prime1));
% plot([0.1 100], [threshold, threshold],'k--')
% plot([0.1 100], [1.36, 1.36],'k--')
% xlim([0.1 100]); ylim([0.1 100]);
% xlabel('Contrast (Michelson)')
% ylabel('Contrast Sensitivity (d prime)'); title('Figure 3')
% set(gca, 'XScale', 'log', 'YScale', 'log')
% set(gca, 'XTick', [0.1, 1, 10, 100], 'XTickLabel', {'0.1', '1', '10', '100'})
% set(gca, 'YTick', [0.1, 1, 10, 100], 'YTickLabel', {'0.1', '1', '10', '100'})


%%

for ii = 1:28, label{ii} = sprintf('%1.3f',expParams.contrastLevels(ii)); end

figure;
hold all;
for ii = 1:length(exParams.eccen)
    plot([0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 6.5, 8],dPrime(end:-1:1,ii));
end
set(gca,'XScale','log');
legend(label);
xlabel('Spatial frequency (cpd)')
ylabel('Sensitivity')

title('CSF eccen 10')
%%
set(gca, 'XTick', [-1, 0, 1, 2], 'XTickLabel', {'0.1', '1', '10', '100'})
set(gca, 'YTick', [-1, 0, 1, 2], 'YTickLabel', {'0.1','1', '10','100'})
box off;
%%

% dPrimeFunction2 = @(alpha, beta) ((nanmean(beta(:)-alpha(:)))./sqrt(nanmean(beta(:))));

% Intensity
% dPrimeFunction3 = @(alpha, beta) (1.36*sqrt(nanmean(beta(:))));

% N = @(alpha, beta) nansum((beta(:)+alpha(:))/2);
% deltaN = @(alpha, beta) nansum(beta(:)-alpha(:));
%
% dPrimeFunction4 = @(alpha, beta) deltaN(alpha, beta)/N(alpha, beta);



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

