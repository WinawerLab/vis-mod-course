function [spacing] = find_spacing(rxy,x,y,spacing_1,spacing_2)
% find the spacing at a particular retinal location [(x,y) coordinate]
spacing = (1/rxy).*sqrt(x.^2*spacing_1.^2+y.^2.*spacing_2);
