function [result,pind] = makeSteerFRs(dims,numLevels,numOrientations,bandwidth)
% makeSteerFRs: makes the frequency responses of the filters for
% a multiscale image transform.
%
% R = makeSteerFRs([M N],numLevels,numOrientations,bandwidth)
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
% DJH 8/96

numbands= numLevels*numOrientations+2;
pind = ones(numbands,1)*dims;

p = numOrientations-1;
Const = sqrt((2^(2*p)*(factorial(p))^2)/(factorial(2*p)*(p+1)));
[f1,f2] = freqspace(dims);
[wx,wy] = meshgrid(f1,f2);
r = sqrt(wx.^2 + wy.^2);
theta = atan2(wy,wx);

result = zeros(numbands*prod(dims),1);

% bandpass bands
for level = 1:numLevels
  for orientation = 0:p;
    thetaOffset = orientation * pi/numOrientations;
    CtrFreq = pi/2^(level*bandwidth);
    band = i^p * Const * cos(theta-thetaOffset).^p ...
	.* logRaisedCos(r,CtrFreq,bandwidth);
    b=(level-1)*numOrientations+(orientation+1)+1;
    indices=pyrBandIndices(pind,b);
    result(indices) = band(:);
  end
end

% hi band
b=1;
hi = logRaisedCosHi(r,pi/(2^bandwidth),bandwidth);
indices=pyrBandIndices(pind,b);
result(indices) = hi(:);

% lo band
b=size(pind,1);
lo = logRaisedCosLo(r,pi/(2^(bandwidth*numLevels)),bandwidth);
indices=pyrBandIndices(pind,b);
result(indices) = lo(:);
