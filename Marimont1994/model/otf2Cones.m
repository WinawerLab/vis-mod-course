% usage:    otf = otf2Cones(params,dispFig)
% by:       Michael Jigo
% purpose:  Generate cone pigment sensitivitiies from the optics.
%
% input:
% p         contains the parameter structure containing the OTF (see inputParams for more info).

function p = otf2Cones(p,dispFig)
if ~exist('dispFig','var')
   dispFig = 0;
end

%% Weight OTF by cone spectral sensitivity
load cones
long = repmat(cones(1,:)',1,size(p.otf,2));
med = repmat(cones(2,:)',1,size(p.otf,2));
short = repmat(cones(3,:)',1,size(p.otf,2));

long = p.otf.*long;
med = p.otf.*med;
short = p.otf.*short;

p.cones(1,:,:) = long;
p.cones(2,:,:) = med;
p.cones(3,:,:) = short;

%% Plot
if dispFig
   figure('Name','Cone sensitivities');
   plots = {long med short};
   plotTitles = 'LMS';
   for i = 1:length(plots)
      subplot(3,2,(i-1)*2+1);
      surf(plots{i}); colormap white
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
      set(gca,'View',[67.3 32.4]);
      title(plotTitles(i));

      subplot(3,2,(i-1)*2+2);
      plot(p.lambda*10^9,plots{i},'k-');
      xlabel('Wavelength (nm)');
      set(gca,'TickDir','out');
      set(gca,'XLim',[min(p.lambda) max(p.lambda)]*10^9);
      title(plotTitles(i));
      box off
   end
end
