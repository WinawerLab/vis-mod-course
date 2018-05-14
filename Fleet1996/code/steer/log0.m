function y = log0(x,base,zero_val)
if ~exist('base')
  base = exp(1);
end
if ~exist('zero_val')
  zero_val = 0;
end
y = zero_val*ones(size(x));
i = find(x > 0);
y(i) = log(x(i)) ./ log(base);

    
