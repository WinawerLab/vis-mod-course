function [grating] = mk_grating(sf,orientation,size_deg)
%create a grating of specified sf/orientation/size in degrees of visual
%angle
deltaX = 0.05; %deg vis angle (1/120)
s = -size_deg:deltaX:size_deg; %space
%%%%
[x,y] = meshgrid(linspace(-4/2,4/2,length(s)+1));
x = x(1:end-1,1:end-1);
y = y(1:end-1,1:end-1);
%make ramp
v_ramp = sin(orientation*pi/180)*x-cos(orientation*pi/180)*y;
%make grating
grating = cos(2*pi*sf*v_ramp);
