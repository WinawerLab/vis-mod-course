function [deg]=mm2deg_ecc(mm)
  
    deg = [3.556*mm+0.05993*mm.^2-0.007358*mm.^3+0.0003027*mm.^4]';
