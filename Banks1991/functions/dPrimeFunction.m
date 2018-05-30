function dprime = dPrimeFunction()

% Define d-prime function
dprime = @(alpha, beta)  (sum( (beta(:)-alpha(:)) .* log(beta(:)./alpha(:))  ) ./ ...
    sqrt(0.5* sum( (beta(:)+alpha(:)) .* log( (beta(:)./alpha(:))).^2 )) );


    
% dPrimeFunction2 = @(alpha, beta) ((nanmean(beta(:)-alpha(:)))./sqrt(nanmean(beta(:))));

% Intensity
% dPrimeFunction3 = @(alpha, beta) (1.36*sqrt(nanmean(beta(:))));

% N = @(alpha, beta) nansum((beta(:)+alpha(:))/2);
% deltaN = @(alpha, beta) nansum(beta(:)-alpha(:));
%
% dPrimeFunction4 = @(alpha, beta) deltaN(alpha, beta)/N(alpha, beta);

