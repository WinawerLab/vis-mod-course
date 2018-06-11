function [sse,fit] = fit_displacementZoneMRGCF(params,y,ecc)
    %this function fits values outside the displacement zone 
    %minimizes sse
    a = params(1);
    r2 = params(2);
    re = params(3);
    drasdo = 1/1.12*(1+(ecc/41.03)).^-1;
    %
    fit = drasdo.*(a*(1+ecc./r2).^-2+(1-a)*exp(-ecc/re));
    sse = sum(y-fit)^2;

