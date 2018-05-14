function [pyr,pind] = buildQuadBands(im, freqRespsImag, freqRespsReal)
% buildQuadBands: builds quadrature pair multiscale subbands,
% given the frequency responses of the quadrature pairs of filters
%
% [pyr,pind] = buildQuadBands(im,freqRespsImag,freqRespsReal)
%    im: input image
%    freqRespsImag and freqRespsReal: filter frequency responses
%       returned by makeQuadFRs.
%
% DJH 8/96

if (size(freqRespsImag)~=size(freqRespsReal))
  error('freqRespsImag and freqRespsImag are incompatible')
end

[pyrImag,pind]=buildSteerBands(im,freqRespsImag);
pyrReal=buildSteerBands(im,freqRespsReal);
pyr=pyrReal+i*pyrImag;

