function [deg_adj] = mm2deg_dens(mm,loc)
    if strcmp(loc,'nasal')
        offset = -1.5; %adjust for offset in visual to optical
    elseif strcmp(loc,'temporal')
        offset =  1.5;
    elseif strcmp(loc,'superior')
        offset = 0.5;
    elseif strcmp(loc,'inferior')
        offset = -0.5
    end
   %convert to mm
   deg = mm2deg_ecc(mm);
   %convert to deg^2
   deg_sqr = 0.0752+(5.846*10^-5)*deg-(1.064*10^-5)*deg.^2+(4.116*10^-8)*deg.^3;
   %adjust for offset
   deg_adj = deg_sqr.*(mm-offset)+(deg_sqr.*offset);
