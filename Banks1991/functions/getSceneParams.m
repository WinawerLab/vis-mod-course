function [sparams, fov2deg] = getSceneParams()


% Scene field of view
sparams.fov       = 1;   % scene field of view in degrees (diameter)
fov2deg           = 1/sparams.fov;
sparams.distance  = 0.57;  % meters

return