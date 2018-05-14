function [cMosaic, emPaths] = getConeMosaic(expParams, thisEccen, deg2m, sparams, ois, whichObserver)

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
%     if whichObserver == 'human'
%         cMosaic.noiseFlag = 'random';
%     else
    cMosaic.noiseFlag = 'none';

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
    
    
    return
