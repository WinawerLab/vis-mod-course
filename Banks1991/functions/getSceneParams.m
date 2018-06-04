function [sparams, fov2deg] = getSceneParams()


% Scene field of view
sparams.fov       = 1;   % scene field of view in degrees (diameter)
fov2deg           = 1/sparams.fov;
sparams.distance  = 0.57;  % meters
sparams.illuminantphotons = 1348 * 10; % retinal illumance, photopic Td. Transformed to photons by multiplying by 10.
sparams.meanluminance = 762; % cd/m2

return