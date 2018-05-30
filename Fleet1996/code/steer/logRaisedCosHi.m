function y = logRaisedCosHi(r,CtrFreq,bandwidth)
CtrFreq = CtrFreq * 2^bandwidth;
Rarg = (pi / bandwidth) * log0(pi/CtrFreq*r,2,-pi);
y = sqrt(1/2*(cos(Rarg)+1));
y = (-pi<Rarg) .* (Rarg<0.0) .* y + (Rarg>=0.0) .* ones(size(y));

  
