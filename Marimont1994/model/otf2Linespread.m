% usage:    otf = otf2Linespread(params,dispFig)
% by:       Michael Jigo
% purpose:  Generate linespread function from Optical Transfer Function.
%
% input:
% p         contains the parameter structure containing the OTF (see inputParams for more info).

function p = otf2Linespread(p,dispFig)
if ~exist('dispFig','var')
   dispFig = 0;
end

%% Compute linespread by transfroming from frequency -> spatial domain
for i = 1:size(p.otf,1)
   otf = p.otf(i,:);
   otf(isnan(otf)) = [];
   p.linespread(i,:) = abs(ifftshift(ifft(otf)));
end
% add nans back to keep matrix size constant
nanIdx = find(all(isnan(p.otf),1));
p.linespread(:,nanIdx) = nan;

%% Plot
if dispFig
   figure('Name','Linespread');
   surf(p.linespread); colormap white
   set(gca,'XTick',linspace(0,size(p.linespread,2),6));
   set(gca,'YTick',linspace(0,p.dLambda,4));
   spaceTicks = 0:0.2:1;
   s = cellfun(@(x) sprintf('%.1f',x),num2cell(spaceTicks),'UniformOutput',0);
   lambdaTicks = round(p.lambda(round(linspace(1,p.dLambda,4)))*1e9*1e3)./1e3;
   l = cellfun(@(x) sprintf('%.0f',x),num2cell(lambdaTicks),'UniformOutput',0);
   set(gca,'XTickLabel',s);
   set(gca,'YTickLabel',l);
   xlabel('degrees');
   ylabel('Wavelength (nm)');
   box off
end
