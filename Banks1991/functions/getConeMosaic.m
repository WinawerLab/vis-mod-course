function cMosaic = getConeMosaic(expParams, thisEccen, deg2m, sparams, ois, whichObserver, segment, verbose)

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
    if whichObserver == 'human'
        cMosaic.noiseFlag = 'random';
    else
        cMosaic.noiseFlag = 'none';
    end

    % Not sure why these have to match, but there is a bug if they don't.
    cMosaic.integrationTime = ois(1).timeStep;

    % There are no eyemovements, but I think you need to have emPaths defined in
    % order to get time varying absorption rates (because it's an oisSequence)
    cMosaic.emGenSequence(ois(1).length, 'nTrials', 1);    
    cMosaic.emPositions = zeros(1,ois(1).length,2);

    % implement th inner segment aperture to correct for proportion covered
    propCovered = getBanks1991ConeCoverage(thisEccen, 'coneSegment', segment, 'verbose', verbose);

    cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered;
    cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered;
    
    
    return
