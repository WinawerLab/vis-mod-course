function sd = table1SizetoSD()

% Reported in table 1
spatialFrequencies  = [.25 .4 .65 1 1.6 2.6 4 6.5 8 10 16 26]; % cycles per degree
eccentricities      = [0 2 5 10 20 40]; % degrees

[sf, eccen] = meshgrid(spatialFrequencies, eccentricities);

numCycles   = [...
    
                NaN .8 1.05 1.2 1.35 1.5 1.75 1.9 NaN 2.03 2.05 2.05; ...   Fovea

                NaN NaN 1 1.3 1.7 2.4 3.2 3.21 NaN 3.79 3.99 4.22; ...      2 deg

                1 1.21 1.51 1.71 2.2 2.5 3.4 4.19 NaN 4.8 NaN NaN; ...      5 deg

                0.8 1 1.49 2 2.39 3.02 3.99 4.5 4.5 NaN NaN NaN; ...        10 deg

                1.7 1.7 1.81 1.90 2 2 2 NaN NaN NaN NaN NaN;  ...           20 deg

                1.7 1.7 1.7 1.7 1.7 NaN NaN NaN NaN NaN NaN NaN; ...        40 deg

                ];

% Calculate the length in degrees of one cyce
cycleLength = 1./sf;          % degrees per cycle

% Calculate how many degrees is 1 SD for each spatial frequency and
% eccentricity
sd = cycleLength .* numCycles/2;