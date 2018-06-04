function [tparams, tseries] = getStimulusParams(thisContrast, targetSize, spatFreq)

% Function to create test and scene params to create Optical Image
% Sequences (OIS) later

% Stimuli are for a 2-AFC detection task. 
% There will always be one mean luminance, gray stimulus and one oriented
% gabor


% Temporal properties of one trial
tStep            = 0.01;                % Time step for optical image sequence (seconds)
tSamples         = (0:tStep:0.20);      % seconds
% timesd           = .01;               % sd of temporal Gaussian window

[~, fov2deg] = getSceneParams();

% Gabor parameters
tparams{1}(1)           = harmonicP;           % Function to get standard Gabor
tparams{1}(1).ang       = pi/2;                % Gabor orientation (radians) - for now, horizontal
tparams{1}(1).freq      = spatFreq*fov2deg;    % Spatial frequency (cycles/FOV)
tparams{1}(1).GaborFlag = targetSize;          % Gaussian window (1 SD)

% Define gabors as test params, and add another one for the blank stimulus
tparams{1}(1).contrast  = 0;                   % Presumably michelson, [0 1]
tparams{1}(2) = tparams{1}(1);
tparams{1}(2).contrast  = 0;

tparams{2} = tparams{1};
tparams{2}(2).contrast = thisContrast;

% Make a Gaussian temporal modulation that brings the stimulus on and off
tseries = zeros(length(tSamples),1);
tseries(5:15) = 1; % HACK
% stimWeights = ieScale(fspecial('gaussian',[1,200],50), 0, 1);