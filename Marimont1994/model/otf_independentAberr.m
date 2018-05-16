% usage:    otf = otf_indAberr(params,dispFig)
% by:       Michael Jigo
% purpose:  Add wavelength-independent aberrations to Optical Transfer Function
%           inputted in "p" structure carried over from "otf.m".

function p = otf_indAberr(p,dispFig)
if ~exist('dispFig','var')
   dispFig = 0;
end

%% Wavelength-independent aberrations scale factor
K = 0.3481+0.6519*exp(-0.1212*p.freq);

%% Add aberrations to OTF
p.otf = p.otf.*K;

%% Plot
if dispFig
   figure('Name',['OTF w/',char(955),'-independent aberrations']);
   surf(p.otf); colormap white
   set(gca,'XTick',linspace(0,length(p.freq),5));
   set(gca,'YTick',linspace(0,p.dLambda,4));
   freqTicks = round(p.freq(round(linspace(1,length(p.freq),5)))*1e3)./1e3;
   f = cellfun(@(x) sprintf('%.0f',x),num2cell(freqTicks),'UniformOutput',0);
   lambdaTicks = round(p.lambda(round(linspace(1,p.dLambda,4)))*1e9*1e3)./1e3;
   l = cellfun(@(x) sprintf('%.0f',x),num2cell(lambdaTicks),'UniformOutput',0);
   set(gca,'XTickLabel',f);
   set(gca,'YTickLabel',l);
   xlabel('Spatial frequency (cpd)');
   ylabel('Wavelength (nm)');
   box off
end
