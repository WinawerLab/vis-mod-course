% usage:    inputParams
% by:       Michael Jigo
% purpose:  Validate any inputted parameters and insert default parameters where empty
%
% input:
% in        Structure with the following fields:
%           dx          number of samples per degree (default: 40) 
%                       NOTE: currently assumes stimulus spans 1 degree
%           
%           dLambda     number of wavelength samples between 400 and 700 nm (default: 31)
%                       NOTE: this variable depends on how finely cone photopigments were sampled in wavelength space
%          
%           q           1x3 element vector containing the variables (q1,q2,q3) for converting defocus to optical path length error (default: [1.7312 0.63346 0.2131]);
%          
%           defocus     dioptric power of the unaccommodated eye (default: 59.9404)
%           
%           pupilRadius radius of the pupil in meters (default: 0.0015)
%
%           degPerM     multiplicative inverse of meters per degree (default: 3434.04)
%
%           im          input 1-dimensional image. currently, the inputted image will be repeated in each color channel (default: step function)
%
%           If no any parameters are missing, the defaults will be used.
%           If any invalid parameters are entered, they will be removed.

function params = inputParams(in)

%% Validate parameters
validParams = {'dx' 'dLambda' 'q'  'defocus' 'pupilRadius' 'degPerM' 'im'};
defaultParams = {40 31 [1.7312 0.63346 0.21310] 59.9404 0.0015 3434.04};
% add in default image (step)
step = ones(1,defaultParams{1});
step(floor(defaultParams{1}/2):end) = 0;
defaultParams{end+1} = step;

switch nargin
   case 0
      % use default parameters
      for p = 1:length(validParams)
         params.(validParams{p}) = defaultParams{p};
      end
   case 1
      % get the inputted parameters and values
      inParams = fieldnames(in);
      for i = 1:length(inParams)
         inVal{i} = in.(inParams{i});
      end	
      % check if inputted parameters are valid, if not remove them and set undefined parameters to default values
      invalid = ~ismember(inParams,validParams);
      if any(invalid)
         warning('The following parameters are invalid and were removed:');
         for i = find(invalid)
            fprintf('%s \n',inParams{i});
         end
         inVal = inVal(~invalid);
         inParams = inParams(~invalid);
      end

      % if any parameters were undefined, use the defaults
      missingParams = ~ismember(validParams,inParams);
      allParams = [inParams; validParams(missingParams)'];
      allVal = [inVal'; defaultParams(missingParams)'];
      for i = 1:length(allParams)
         params.(allParams{i}) = allVal{i};
      end
end

%% Finalize params ouput
params.space = -0.5:1/params.dx:(0.5-1/params.dx);
params.lambda = linspace(400,700,params.dLambda)*1e-9;
params.freq = (0:params.dx-1)*params.dx/params.dx;
