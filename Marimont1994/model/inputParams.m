% usage:    inputParams
% by:       Michael Jigo
% purpose:  Validate any inputted parameters and insert default parameters where empty

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
