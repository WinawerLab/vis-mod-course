function [freqRespsImag,freqRespsReal,pind] = ...
    makeQuadFRs(dims,numLevels,numOrientations,bandwidth)
% makeQuad: quadrature pairs of makeSteerFRs.
%
% R = makeQuadFRs([M N],numLevels,numOrientations,bandwidth)
%    [M N]: image size
%    numLevels: number of levels/scales.
%    numOrientations: number of orientation subbands at each scale.
%    bandwidth: spatial frequency bandwidth in octaves
%
% Return values
%    result: vector that contains all of the subband images.
%    pind: contains the sizes of each of the resulting subband
%       images, following conventions used by Eero Simoncelli's
%       pyramid code. 
%
% DJH, 8/96

if isOdd(numOrientations)
  error('numOrientations must be even in makeQuadFRs.');
end

numbands= numLevels*numOrientations+2;
p = numOrientations-1;

[freqRespsImag,pind]= ...
    makeSteerFRs(dims,numLevels,numOrientations,bandwidth);

freqRespsReal=freqRespsImag;
for b=2:(numbands-1)
  indices=pyrBandIndices(pind,b);
  freqRespsReal(indices)=abs(freqRespsReal(indices));
end

% Finally, zero out low and high band of imaginary FRs
b=1;
indices=pyrBandIndices(pind,b);
freqRespsImag(indices)=zeros(size(indices));
b=numbands;
indices=pyrBandIndices(pind,b);
freqRespsImag(indices)=zeros(size(indices));


