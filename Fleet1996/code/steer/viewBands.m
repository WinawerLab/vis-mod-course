function viewBands(pyr, pind, wait, rng)
% viewBands: view the subbands of the transform
%
% viewBands( pyr, pind, wait)
%    pyr: the image transform
%    pind: returned by buildQuadBands or buildSteerBands,
%       specifies the sizes of the subbands, following
%       conventions used by Eero Simoncelli's pyramid code.
%    wait: time to wait between displaying successive subband images.
%
% DJH, 8/96

if ~exist('wait')
  wait=1;
end

if ~exist('rng')
  rng = 'auto1';
end

M=pind(1,1);
N=pind(1,2);
totalsize = M*N;
numbands = size(pyr,1)/totalsize;

for b=1:numbands
  displayImage(reshape(pyr(pyrBandIndices(pind,b)),pind(b,1),pind(b,2)),rng);
  pause(wait);
end

