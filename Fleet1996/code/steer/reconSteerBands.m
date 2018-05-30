function result = reconSteerBands( pyr, freqResps, pind )
% reconSteerBands: reconstruct an image from the subband transform.
%
% result = reconSteerBands(pyr,freqResps,pind)
%    pyr: the image transform
%    freqResps: filter frequency responses returned by makeSteerFRs.
%    pind: returned by buildQuadBands or buildSteerBands,
%       specifies the sizes of the subbands, following
%       conventions used by Eero Simoncelli's pyramid code.

% ToDo:
% - check that numbands is integer 
% - check that size(pyr)==size(freqResps)

M=pind(1,1);
N=pind(1,2);
totalsize = M*N;
numbands = size(freqResps,1)/totalsize;

result = zeros(M,N);
for b=1:numbands
  indices=pyrBandIndices(pind,b);
  freqResp = reshape(freqResps(indices),pind(b,1),pind(b,2));
  band = reshape(pyr(indices),pind(b,1),pind(b,2));
  freqBand = fftshift(fft2(band));
  result = result + real(ifft2(fftshift(freqBand.*conj(freqResp))));
end


