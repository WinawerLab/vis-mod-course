function [sse,fit] = fit_displacementZone(params,y,ecc)
    %this function fits values outside the displacement zone 
    %minimizes sse
    dgf = 33163.2; %taken from paper
    a = params(1);
    r2 = params(2);
    re = params(3);
    %
    fit = dgf*(a*(1+ecc./r2).^-2+(1-a)*exp(-ecc/re));
    sse = sum(y-fit)^2;

