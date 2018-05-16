% usage:    otf = otf(params,dispFig)
% by:       Michael Jigo
% purpose:  Create the Optical Transfer function (OTF) given the input 
%           parameters (p).

function p = otf(p,dispFig)
if ~exist('dispFig','var')
   dispFig = 0;
end

%% Spatial frequency (SF) --> Reduced SF
for l = 1:length(p.lambda)
   reducedSF(l,:) = p.degPerM*(p.lambda(l)/(p.defocus*p.pupilRadius)).*p.freq;
end

%% Defocus (D) --> optical path length error (w20)
D = p.q(1)-p.q(2)./((p.lambda*1e6)-p.q(3));
w20 = (p.pupilRadius^2)/2*((p.defocus*D)./(p.defocus+D));

%% Create OTF
% initialize otf field
p.otf = nan(length(p.lambda),size(reducedSF,2));
for l = 1:length(p.lambda)
   for m = 1:size(reducedSF,2)
      intgCount = 1;
      a = (4*pi)/p.lambda(l)*w20(l)*abs(reducedSF(l,m));

      % integrate
      step = 1/p.dx; lim = sqrt(1-(reducedSF(l,m)/2)^2); % limit of integration
      intg = nan(1,length(1/step:1/step:lim));
      for i = step:step:lim
         intg(intgCount) = sin(a*(sqrt(1-i^2)-abs(reducedSF(l,m))/2));
         intgCount = intgCount+1;
      end

      % compute otf
      p.otf(l,m) = 4/(pi*a)*(nansum(intg)/length(intg));
      clear intg
   end
end

%% Plot
if dispFig
   figure('Name','OTF');
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
