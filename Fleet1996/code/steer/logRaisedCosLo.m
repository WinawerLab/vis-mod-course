function y = logRaisedCosLo(r,CtrFreq,bandwidth)
CtrFreq = CtrFreq/2^bandwidth;
Rarg = (pi / bandwidth) * log0(pi/CtrFreq*r,2,0);
y = sqrt(1/2*(cos(Rarg)+1));
y = (0.0<Rarg).*(Rarg<pi).*y + (Rarg<=0.0).*ones(size(y));

  
