function y = logRaisedCos(r,CtrFreq,bandwidth)
Rarg = (pi / bandwidth) * log0(pi/CtrFreq*r,2,bandwidth);
y = sqrt(1/2*(cos(Rarg)+1));
y = (-pi<Rarg) .* (Rarg<pi) .* y;

  
