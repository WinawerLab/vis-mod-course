function disparityMap = disparitySim(I_L,I_R,filterBank,plotMode)
%   Title:   Implementation and expansion of Fleet, Wagner, and Heeger, Vision Research (1995),
%            "Neural Encoding of Binocular Disparity: Energy Models, Position Shifts
%            and Phase Shifts"
%   Purpose: This code expands the model to two dimensions and
%            takes real stereo pairs as input (vs. simple contrived
%            stimuli). Created for Jon Winawer's Computational Modelling of
%            Visual Circuits Course at NYU, Dept. of Psychology, 2018.
%   Author:  Charlie Burlingham
%   Date:    April 21, 2018
%   Depend.: steerable pyramid toolbox, image processing toolbox.
%   Inputs:
%           I_L: left eye image (grayscale)
%           I_R: right eye image(grayscale)
%           filterBank: 0 = gabors, 1 = steerable pyramid
%           plotMode: 0 = don't plot figs, 1 = plot figs
%   Outputs:
%           A disparity map d(x,y)
%
%   Input image examples
% load in real stereo pair
%{
    I_R = imread('im0.png'); % right eye image
    I_R = rgb2gray(I_R);
    I_L = imread('im1.png'); % left eye image is same as right eye image
    I_L = rgb2gray(I_L);

    % load in artifical rings image
    I_L = imread('concentric-circles-paul-sober.jpg'); % right eye image
    disparityPixels = 20;
    I_R = circshift(I_L,disparityPixels,2); % circshift horizontally by ground truth disparity amount
    
    % or use circshifted uniform white noise images
    sizeI = 500;
    I_L = rand(sizeI,sizeI);
    disparityPixels = 2;
    I_R = circshift(I_L,disparityPixels,2); % circshift horizontally by ground truth disparity amount

    % Small Box in Depth 
    sizeI = 500;
    sizeBox = 200;
    disparityMapGT = -5*ones(sizeI,sizeI);
    disparityMapGT(150:350-1,150:350-1) =  10*ones(sizeBox,sizeBox); 
    [I_L,I_R] = stereogram(disparityMapGT,[],0); % don't plot

    subplot(1,2,1); imagesc(I_L); colormap gray; subplot(1,2,2); imagesc(I_R);colormap gray;

    % Gradient 
    sizeI = 500;
    sizeBox = 200;
    disparityMapGT = -5*ones(sizeI,sizeI);
    disparityMapGT(150:350-1,150:350-1) = repmat(linspace(-14,14,sizeBox),[sizeBox 1]); % gradient
    [I_L,I_R] = stereogram(disparityMapGT,[],0); % don't plot
    subplot(1,2,1); imagesc(I_L); colormap gray; subplot(1,2,2); imagesc(I_R);colormap gray;





% figure out how to encode the disparity map explicitly in the scale
factors on the input to stereogram.m ....



    % idk how to generate the stereograms myself, figure this out...
    %{
    I_L = rand(sizeI,sizeI);
    I_R = circshift(I_L,-2,2); % circshift horizontally by ground truth disparity amount

    box = rand(sizeBox,sizeBox);
    
    disparityPixels = 10;
    I_L(150:350-1,150-disparityPixels:350-disparityPixels-1) = box*.9;

    I_R(150:350-1,150+disparityPixels:350+disparityPixels-1) = box*.9;
    %}

%}

%% Setup defaults & clean up working space
addpath(genpath('steer'))
close all; clc;
if ieNotDefined('filterBank'); filterBank = 1; end % steerable pyr default
if ieNotDefined('plotMode'); plotMode = 0; end % don't plot by default

%% Get "ground truth" disparity map, for visualization

disparityRange = [-6 10];
disparityMapGTRecovered = disparity(I_L, I_R, 'BlockSize', 15, ...
    'DisparityRange', disparityRange);

%% plot input images
if plotMode
    figure(1)
    subplot(1,4,1)
    imagesc(I_L)
    axis square
    colormap gray
    subplot(1,4,2)
    imagesc(I_R)
    colormap gray
    axis square
    subplot(1,4,3)
    imagesc(disparityMapGT)
    colormap gray
    axis square
    subplot(1,4,4)
    imagesc(disparityMapGTRecovered)
    colormap gray
    axis square
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
        subplot(1,3,1)
        for kk = 1:size(E_x,3); imagesc(E_x(:,:,kk)); colormap gray; pause(.05);end % uncomment to plot energy for stereo pair\\
        subplot(1,3,1)
        imagesc(pooled_E_x); colormap gray;
        title('Zero Disparity Model','Fontsize',20)
        
    elseif filterBank == 1 % steerable pyramids
        % Look at the filters:
        figure(2)
        displayImage(accessSteerBand(freqRespsReal,pind,numOrientations,1,1));
        viewBands(freqRespsReal,pind,1/4,'auto1');
        
        figure(3)
        % look at the pyramid responses to the image (just to left eye image)
        for ii = 1:size(amp_L,3)
            imagesc(amp_L(:,:,ii)); pause(.1);  %  plot magnitude of quad resps
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Position shift model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = disparityMapGT; % range of possible disparities. the amount the eye-specific input images are actually shifted by in pixels.
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
disparityMap_posShift = s(maxDisparityInd_posShift);

E_x_eq5 = amp_L.^2 + amp_R.^2  + 2.*amp_L.*amp_R .* cos(phase_L-phase_R); % eq 5, this is a version just using the responses to the two shifted images... not assuming they are merely translated copies

if plotMode
    figure(3);
    subplot(1,3,2)
    for ii = 1:size(pooled_E_x_posShift,3); imagesc(pooled_E_x_posShift(:,:,ii)); colormap gray; pause(.1);end % uncomment to plot energy for stereo pair
    
    figure
    imagesc(disparityMap_posShift); colormap gray;
    title('Position Shift Model','Fontsize',20)
end

figure

imagesc(sum(E_x_eq5,3))

 mse(disparityMap_posShift,disparityMapGT)
 mse(disparityMapGT,disparityMapGTRecovered)

% Phase shift model DOESN'T WORK AT ALL YET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase shift model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaPhase = -pi:pi/50:pi;

pooled_E_x_phaseShift = NaN(size(I_L,1),size(I_L,2),length(deltaPhase));
% now sum the power from the two eyes' filter banks
for ii = 1:length(deltaPhase)
    current_deltaPhase = repmat(deltaPhase(ii), size(d2));
    E_x_phaseShift = 2.*amp_L.^2 + ( 2.*amp_L.^2 .* cos((d2.*k) - current_deltaPhase) )  ; % eq 10, assuming the inputs are just horizontally translated copies of one another
    
    pooled_E_x_phaseShift(:,:,ii) = sum(E_x_phaseShift,3); % sum across orientation and spatial frequency channels
end


[maxEnergy_phaseShift maxDisparityInd_phaseShift] = max(pooled_E_x_phaseShift,[],3); % find max along the third dimension, representing a "max" readout of the Energy as a function of disparity at each pixel.
disparityMap_phaseShift = deltaPhase(maxDisparityInd_phaseShift);


if plotMode
    figure(3);
    for ii = 1:size(pooled_E_x_phaseShift,3); imagesc(pooled_E_x_phaseShift(:,:,ii)); colormap gray; pause(.1);end % uncomment to plot energy for stereo pair
    
    figure
    imagesc(disparityMap_phaseShift); colormap gray;
    title('Phase Shift Model','Fontsize',20)
end



end