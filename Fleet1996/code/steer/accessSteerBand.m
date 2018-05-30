function res = accessSteerBand(pyr,pind,numOrientations,level,orientation);
% accessSteerBand: pulls specified band (level,orientation) out
% of a subband transform "data structure".
%
% res = accessSteerBand(pyr,pind,numOrientations,level,orientation);
%    pyr: the image transform
%    pind: returned by buildQuadBands or buildSteerBands,
%       specifies the sizes of the subbands, following
%       conventions used by Eero Simoncelli's pyramid code.
%    numOrientations: number of orientation subbands at each scale.
%    level: desired level/scale.
%    orientation: desired orientation.
%    
% DJH 8/96

b=(level-1)*numOrientations+orientation+1;
res=reshape(pyr(pyrBandIndices(pind,b)),pind(b,1),pind(b,2));
