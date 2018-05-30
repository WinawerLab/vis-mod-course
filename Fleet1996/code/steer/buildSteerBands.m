function [pyr,pind] = buildSteerBands(im, freqResps)
% buildSteerBands: builds subbands multiscale of a multiscale
% image transform, given the frequency responses of the filters. 
%
% [pyr,pind] = buildSteerBands(im,freqResps)
%    im: input image
%    freqResps: filter frequency responses returned by makeSteerFRs.
%
% DJH 8/96


% ToDo:
% - check that numbands is integer 
% - check that size(pyr)==size(freqResps)

[M N] = size(im);
totalsize = M*N;
numbands = size(freqResps,1)/totalsize;
pind = ones(numbands,1)*size(im);

fourier = fftshift(fft2(im));

pyr = zeros(size(freqResps));
for b=1:numbands
  indices=pyrBandIndices(pind,b);
  freqResp = reshape(freqResps(indices),pind(b,1),pind(b,2));
  band = real(ifft2(fftshift(fourier.*freqResp))); % CSB: convolution
  pyr(indices) = band(:);
end


