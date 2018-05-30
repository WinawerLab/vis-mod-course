function lev = maxLevel(dims,bandwidth)
% lev = maxLevel(dims,bandwidth)
% DJH 8/96

lev = floor((log(min(dims))/log(2)-2)/bandwidth);
