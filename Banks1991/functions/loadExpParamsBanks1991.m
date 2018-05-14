function expParams = loadExpParamsBanks1991()


% Get std's of Gaussian window
[sd, sf, eccen, cycleLength] = table1SizetoSD;

expParams.sd          = sd;
expParams.sf          = sf;
expParams.eccen       = eccen;
expParams.cycleLength = cycleLength;

% General parameters
expParams.nTrials           = 1;                            % Number of trials per stimulus condition
expParams.contrastLevels    = [0.001:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1];   % Contrast log levels
expParams.polarangle        = 0;                            % Horizontal meridian
expParams.eyemovements      = [0 0 0]';                     % No eye movments
expParams.verbose           = true;

% cMosaic Params
expParams.whichEye           = 'left';

return