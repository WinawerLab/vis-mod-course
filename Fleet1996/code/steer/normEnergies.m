function result = normEnergies(pyr,pind,numOrientations,sigma)
% normEnergies: computes contrast normalized local energy
% responses from a quadrature subband transform.
%
% result = normEnergies(pyr,pind,numOrientations,sigma)
%    pyr: the image transform
%    pind: returned by buildQuadBands or buildSteerBands,
%       specifies the sizes of the subbands, following
%       conventions used by Eero Simoncelli's pyramid code.
%    numOrientations: number of orientation subbands at each scale.
%    sigma: semi-saturation constant for the contrast normalization.
%
% DJH, 8/96

if ~exist('sigma')
  sigma=0;
end

dims=pind(1,:);
numLevels =(size(pind,1)-2)/numOrientations;

result=zeros(size(pyr));

energies = real(pyr).^2 + imag(pyr).^2;
normalizers = zeros((numLevels+1)*prod(dims),1);
normPind = ones(numLevels+1,1)*dims;

% Compute hipass normalizer (spatial blurring here because there
% is no quadrature pair).
hi=blur(pyrHi(energies,pind));
normIndices=pyrBandIndices(normPind,1);
normalizers(normIndices)=hi(:);

% Compute normalizers
for level = 1:numLevels
  for orientation = 1:numOrientations
    b=(level-1)*numOrientations+orientation+1;
    pyrIndices=pyrBandIndices(pind,b);
    normIndices=pyrBandIndices(normPind,level+1);
    normalizers(normIndices)=normalizers(normIndices)+energies(pyrIndices);
  end
end

% Normalize each bandpass energy band, except skip the lowest sf.
for level = 1:(numLevels-1)
  for orientation = 1:numOrientations
    b=(level-1)*numOrientations+orientation+1;
    pyrIndices=pyrBandIndices(pind,b);
    normIndices1=pyrBandIndices(normPind,level);
    normIndices2=pyrBandIndices(normPind,level+1);
    normIndices3=pyrBandIndices(normPind,level+2);
    result(pyrIndices)=energies(pyrIndices) ./ ...
	(normalizers(normIndices1)+normalizers(normIndices2) + ...
	normalizers(normIndices3)+sigma^2);
  end
end

