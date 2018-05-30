function [errorD disp] = disparitySim(I_L,I_R,filterBank,plotMode)
%   Title:   Implementation and expansion of Fleet, Wagner, and Heeger, Vision Research (1995),
%            "Neural Encoding of Binocular Disparity: Energy Models, Position Shifts
%            and Phase Shifts"
%   Purpose: This code expands the model to two dimensions and
%            takes real stereo pairs as input (vs. simple contrived
%            stimuli). Created for Jon Winawer's Computational Modelling of
%            Visual Circuits Course at NYU, Dept. of Psychology, 2018.
%   Author:  Charlie Burlingham
%   Date:    May 14, 2018
%   Dependencies: steerable pyramid toolbox, image processing toolbox.
%   Inputs:
%           I_L: left eye image (grayscale)
%           I_R: right eye image(grayscale)
%           disp.disparityMapGT: A ground truth disparity map d(x,y)
%           filterBank: 0 = gabors, 1 = steerable pyramid
%           plotMode: 0 = don't plot figs, 1 = plot figs
%           
%   Outputs:
%           disp: struct with GT disparity map d(x,y) and decoded/recovered
%           versions d_hat(x,y) from each model
%           errorD: struct with mse between ground truth disparity map
%           and each model estimate.
%
%   Notes: Gabor filter bank not supported fully yet. Use steerable pyramid instead!
%           
%   Example Model Inputs:
%{


    % Recommend using these:

    % Small Box in Depth 
    sizeI = 500;
    sizeBox = 200;
    disp.disparityMapGT = -5*ones(sizeI,sizeI);
    disp.disparityMapGT(150:350-1,150:350-1) =  10*ones(sizeBox,sizeBox); 
    [I_L,I_R] = stereogram(disp.disparityMapGT,[],0); % don't plot
    subplot(1,2,1); imagesc(I_L); colormap gray; subplot(1,2,2); imagesc(I_R);colormap gray;

    % Gradient box
    sizeI = 500;
    sizeBox = 200;
    disp.disparityMapGT = -3*ones(sizeI,sizeI);
    disp.disparityMapGT(150:350-1,150:350-1) = repmat(linspace(-5,5,sizeBox),[sizeBox 1]); % gradient
    [I_L,I_R] = stereogram(disp.disparityMapGT,[],0); % don't plot
    subplot(1,2,1); imagesc(I_L); colormap gray; subplot(1,2,2); imagesc(I_R);colormap gray;


    % Still in testing:

    I_R = imread('view1.png'); % right eye image
    I_R = rgb2gray(I_R);
    I_L = imread('view5.png'); % left eye image is same as right eye image
    I_L = rgb2gray(I_L);

    disp.disparityMapGT = double(imread('disp1.png')); % load GT disparity map,
    % ATTN: this is a problem because this is only the one for the right
    %eye image... anyway, for now going to assume it's the canonical map.

    I_L = I_L(1:350,1:350); I_R = I_R(1:350,1:350); % crop
    disp.disparityMapGT = disp.disparityMapGT(1:350,1:350); 
    % ATTN: only crop for the middlebury 2006 benchmark.... the x and y 
    %lengths need to be the same for the input image 

    % load in artifical rings image
    I_L = imread('concentric-circles-paul-sober.jpg'); % right eye image
    disparityPixels = 20;
    I_R = circshift(I_L,disparityPixels,2); % circshift horizontally by ground truth disparity amount
    
    % or use circshifted uniform white noise images
    sizeI = 500;
    I_L = rand(sizeI,sizeI);
    disparityPixels = 2;
    I_R = circshift(I_L,disparityPixels,2); % circshift horizontally by ground truth disparity amount

%}

%% Setup defaults & clean up working space
addpath(genpath('steer'))
close all; clc;
if ieNotDefined('filterBank'); filterBank = 1; end % steerable pyr default
if ieNotDefined('plotMode'); plotMode = 1; end % plot by default

%% Get "ground truth" disparity map, for visualization

disparityRange = [-6 10]; % default
%disparityRange = [-2 158]; % use this for natural images whose d(x,y) uses pixel intensities as proxy for disparity 


disp.disparityMapGTRecovered = disparity(I_L, I_R, 'Method', 'SemiGlobal','BlockSize', 19, ...
    'DisparityRange', disparityRange);

%% plot input images
if plotMode
    figure(1)
    subplot(1,4,1)
    imagesc(I_L)
    axis square
    colormap gray
    title('Left Eye Input','Fontsize',20)
    subplot(1,4,2)
    imagesc(I_R)
    colormap gray
    axis square
    title('Right Eye Input','Fontsize',20)
    subplot(1,4,3)
    imagesc(disp.disparityMapGT)
    colormap gray
    axis square
    title('Ground Truth Disparity Map D(x,y)','Fontsize',20)
    subplot(1,4,4)
    imagesc(disp.disparityMapGTRecovered)
    colormap gray
    axis square
    title('SemiGlobal Recovered D_hat(x,y)','Fontsize',20)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero disparity model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numOrientations = 10;
numSF = 9;
SFBandwidth = 1;
if filterBank == 0
    % create filter bank
    lambda = linspace(2,38,numSF); % wavelength; controls spatial frequency and also size jointly
    ori = linspace(0,160,numOrientations);
    AR = 1; % isometric Gabor
    filt_L = gabor(lambda,ori,'SpatialFrequencyBandwidth',SFBandwidth,'SpatialAspectRatio',AR); % create filter bank with all combos of orientation and SF
    filt_R = filt_L; % zero disparity model assumes the filters in left and right eye are the same!
    % apply filter bank to eye-specific images separately
    [rho_L, phase_L] = imgaborfilt(I_L,filt_L);
    [rho_R, phase_R] = imgaborfilt(I_R,filt_R);
    
    % now sum the power from the two eyes' filter banks
    E_x = 2.*(rho_L).^2 + 2.*(rho_L).^2 ; % eq. 5. Gives the binocular energy in terms of the amplitude and phase of the monocular response, with the full squaring non-linearity. cos(phase_l - phase_r) = 1 so that term vanishes
    pooled_E_x = sum(E_x,3); % sum across orientation and spatial frequency channels
    
elseif filterBank == 1
    dims = size(I_L); % ATTN: don't hard code this, use the max SF of the gabor filter bank
    numLevels = numSF;
    [freqRespsImag, freqRespsReal, pind]=makeQuadFRs(dims,numLevels,numOrientations,SFBandwidth);
    
    [pyr_L,pind_L] = buildQuadBands(I_L, freqRespsImag, freqRespsReal);
    [pyr_R,pind_R] = buildQuadBands(I_R, freqRespsImag, freqRespsReal);
    
    n = 1;
    for ii = 1:numOrientations
        for jj = 1:numLevels
            res_L(:,:,n) = accessSteerBand(pyr_L,pind_L,numOrientations,jj,ii); % complex
            res_R(:,:,n) = accessSteerBand(pyr_R,pind_R,numOrientations,jj,ii);
            n = n + 1;
        end
    end
    
    % get amplitude (rho) and phase (phi) of complex valued filter responses
    amp_L = abs(res_L); amp_R = abs(res_R);
    phase_L = angle(res_L); phase_R = angle(res_R);
    
    % now sum the power from the two eyes' filter banks
    E_x = 2.*(amp_L).^2 + 2.*(amp_L).^2 ; % eq. 5. Gives the binocular energy in terms of the amplitude and phase of the monocular response, with the full squaring non-linearity. cos(phase_l - phase_r) = 1 so that term turns into another 1 * 2*rho^2 !
    pooled_E_x = sum(E_x,3); % sum across orientation and spatial frequency channels
    
end

if plotMode
    if filterBank == 0 % gabors
        % plot real component of filters in filter bank
        figure(2)
        for ii = 1:length(filt_L)
            if mod(ii,2) == 0 % only plot a subset of filter
                subplot(length(ori),length(lambda),ii)
                imagesc(real(filt_L(1,ii).SpatialKernel)); colormap gray
            end
        end
        suptitle('Real component of Gabor filters as function of SF (x-axis) & orientation (y-axis)')
        figure(3)
        for kk = 1:size(E_x,3); imagesc(E_x(:,:,kk)); colormap gray; pause(.05);end % uncomment to plot energy for stereo pair\\
        subplot(1,3,1)
        imagesc(pooled_E_x); colormap gray;
        title('Zero Disparity Model, Energy E(x,y) Pooled Across SF and Ori','Fontsize',20)
        
    elseif filterBank == 1 % steerable pyramids
        % Look at the filters:
        figure(2)
        displayImage(accessSteerBand(freqRespsReal,pind,numOrientations,1,1));
        viewBands(freqRespsReal,pind,1/4,'auto1');
        
        figure(3)
        % look at the pyramid responses to the image (just to left eye image)
        for ii = 1:size(amp_L,3)
            imagesc(amp_L(:,:,ii)); pause(.1); colormap gray;  %  plot magnitude of quad resps
        end
        
        figure 
        imagesc(pooled_E_x); colormap gray; title('Zero Disparity Model, Energy E(x,y) Pooled Across SF and Ori','Fontsize',20)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Position shift model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = disp.disparityMapGT; % range of possible disparities. the amount the eye-specific input images are actually shifted by in pixels.
d2 = repmat(d,[ 1 1 size(phase_L,3) ]); % just copied along third dimension, so we can multiply with k.
s = -15:.5:15; % shift amount: the amount in pixels the centers of the eye-specific filters/RFs are shifted by. Range should be constrained about what we know about the disparities present in input images.

% compute  gradient of response phase "image" to get instantaneous frequency image
for ii = 1:size(phase_L,3)
    [k(:,:,ii), ~] = imgradient(phase_L(:,:,ii),'prewitt'); % k is instantaneous frequency.       ATTN: is prewitt the best thing to use here? the gradients look weird...
end

pooled_E_x_posShift = NaN(size(I_L,1),size(I_L,2),length(s));
% now sum the power from the two eyes' filter banks
for ii = 1:length(s)
    current_s = repmat(s(ii), size(d2));
    E_x_posShift = 2.*amp_L.^2 + ( 2.*amp_L.^2 .* cos(k.*(d2 - current_s)) ) ; % eq 10, assuming the inputs are just horizontally translated copies of one another
    
    pooled_E_x_posShift(:,:,ii) = sum(E_x_posShift,3); % sum across orientation and spatial frequency channels
end

[maxEnergy_posShift maxDisparityInd_posShift] = max(pooled_E_x_posShift,[],3); % find max along the third dimension, representing a "max" readout of the Energy as a function of disparity at each pixel.
disp.disparityMap_posShift = s(maxDisparityInd_posShift);

E_x_eq5 = amp_L.^2 + amp_R.^2  + 2.*amp_L.*amp_R .* cos(phase_L-phase_R); % eq 5, this is a version just using the responses to the two shifted images... not assuming they are merely translated copies

if plotMode
    figure(3);
    for ii = 1:size(pooled_E_x_posShift,3); imagesc(pooled_E_x_posShift(:,:,ii)); colormap gray; pause(.1);end % uncomment to plot energy for stereo pair
    
    figure
    plot(s,squeeze(pooled_E_x_posShift(20,30,:))); xlabel('Position Shifts (Pixels)'); ylabel('Energy E(20,30) Pooled over SF and Ori'); title('Example Position Shift Disparity Preference','Fontsize',20)
    
    figure
    imagesc(disp.disparityMap_posShift); colormap gray;
    title('Position Shift Model Decoded Disparity Map D_{hat}(x,y)','Fontsize',20)
    
end



%  figure; imagesc(sum(E_x_eq5,3));

errorD.msePosShift = mse(disp.disparityMap_posShift,disp.disparityMapGT) % error between GT d(x,y) and decoded/recovered d_hat(x,y) according to position shift model
 
errorD.mseGTdefaultAlg =  mse(disp.disparityMapGT,disp.disparityMapGTRecovered) % error between GT d(x,y) and decoded/recovered d_hat(x,y) according to some standard matlab algorithm
 
 % In future could compare more models, although it's not fair because
 % Fleet model has explicit access to the GT disparity map and the others
 % do not. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase shift model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaPhase = -3*pi:pi/20:3*pi;

pooled_E_x_phaseShift = NaN(size(I_L,1),size(I_L,2),length(deltaPhase));
% now sum the power from the two eyes' filter banks
for ii = 1:length(deltaPhase)
    current_deltaPhase = repmat(deltaPhase(ii), size(d2));
    E_x_phaseShift = 2.*amp_L.^2 + ( 2.*amp_L.^2 .* cos((d2.*k) - current_deltaPhase) )  ; % eq 10, assuming the inputs are just horizontally translated copies of one another
    
    pooled_E_x_phaseShift(:,:,ii) = sum(E_x_phaseShift,3); % sum across orientation and spatial frequency channels
end


[maxEnergy_phaseShift maxDisparityInd_phaseShift] = max(pooled_E_x_phaseShift,[],3); % find max along the third dimension, representing a "max" readout of the Energy as a function of disparity at each pixel.
disp.disparityMap_phaseShift = deltaPhase(maxDisparityInd_phaseShift);


if plotMode
    figure(3);
    for ii = 1:size(pooled_E_x_phaseShift,3); imagesc(pooled_E_x_phaseShift(:,:,ii)); colormap gray; pause(.1); end % uncomment to plot energy for stereo pair
    
    figure; 
    plot(deltaPhase,squeeze(pooled_E_x_phaseShift(1,1,:))); ylabel('Energy E(1,1)'); xlabel('Phase Shift Amount (Radians)') % Energy at pixel 1,1 across disparities shows that the solution is underconstrained. there are peaks periodically corresponding to phase wrap around
    title('solution for disparity pref of phase shift neurons: underconstained')
    
    figure
    imagesc(disp.disparityMap_phaseShift); colormap gray;
    title('Phase Shift Model Decoded Disparity Map D_{hat}(x,y)','Fontsize',20)
end

errorD.msePhaseShift = mse(disp.disparityMap_phaseShift,disp.disparityMapGT) % error between GT d(x,y) and decoded/recovered d_hat(x,y) according to phase shift model
 
errorD.mseGTdefaultAlg =  mse(disp.disparityMapGT,disp.disparityMapGTRecovered) % error between GT d(x,y) and decoded/recovered d_hat(x,y) according to some standard matlab algorithm
 

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hybrid Position and Phase Shift Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The encoding model is straightforward, but decoding without ambiguity
% will require pooling only the Energy responses only from neurons with
% specific SF preferences... 

% To Do: Can do  SF dependent pooling of Energy for the phase shift model
% above and that should help improve the decoding considerably


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Disparity Tuning Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This will require contriving stimuli with different GT disparities and
% simulating the responses of individual energy neurons... Should be very
% simple, just haven't had time to code it up yet. Will illustrate the
% false peaks phenomenon described in the Fleet et al. 1996 paper.



end