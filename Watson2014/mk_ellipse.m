function [x,y] = mk_ellipse(rx,ry,x0,y0)
% rx = horizontal radius
% ry = vertical radius
% x0 = ellipse center x0
% y0 = ellipse center y0
t = -pi:0.01:pi;
x = x0+rx*cos(t); y = y0+ry*sin(t);