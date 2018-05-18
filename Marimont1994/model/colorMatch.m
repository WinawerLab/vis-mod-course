% usage:    colorMatch(im,params,dispFig)
% by:       Michael Jigo
% purpose:  Generate matching color images after taking chromatic aberration
%           into account. Currently works only for 1D images.

function colorMatch(p,disp1,disp2,dispFig)
if ~exist('dispFig','var')
   dispFig = 0;
end

%% Get Fourier image
sfs = p.freq(1:floor(length(p.freq)/2));
if size(p.im,1)==1
   p.im = repmat(p.im,3,1);
end
for i = 1:size(p.im,1)
   freqIm(i,:) = fft(p.im(i,:),p.dx)*(1/p.dx);
end

%% Compute calibration matrix
C1 = nan(length(sfs),size(p.cones,1),size(disp1,1));
C2 = C1;
for s = 1:length(sfs)
   for c = 1:size(p.cones,1)
      for pr = 1:size(disp1,1)
         % disp1 --> cones
         C1(s,c,pr) = nansum(disp1(pr,:).*p.cones(c,:,s));
         % disp2 --> cones
         C2(s,c,pr) = nansum(disp2(pr,:).*p.cones(c,:,s));
      end
   end
end

%% Photopigment resposnes (Display 1)
pigFreq = nan(size(p.cones,1),length(sfs));
for s = 1:length(sfs)
   % extract amplitude at frequency s
   tempIm = freqIm(:,s);

   % photopigment response (frequency domain)
   pigFreq(:,s) = squeeze(C1(s,:,:))*tempIm;
end
% photopigment response (spatial domain)
for c = 1:size(p.cones,1)
   pigSpace(c,:) = real(ifft(pigFreq(c,:),p.dx));
end

%% Predict photopigment responses (Display 2)
pigFreq2 = nan(size(p.cones,1),length(sfs));
for s = 1:length(sfs)
   pigFreq2(:,s) = pinv(squeeze(C2(s,:,:)))*pigFreq(:,s);
end
for c = 1:size(p.cones,1)
   pigSpace2(c,:) = real(ifft(pigFreq2(c,:),p.dx));
end

if dispFig
   col = 'rgb';
   %% Display phosphors
   %figure('Name','Display primaries');
   %disps = {disp1 disp2};
   %for d = 1:length(disps)
   %   subplot(1,2,d);
   %   for j = 1:size(disps{d},1)
   %      plot(p.lambda*10^9,disps{d}(j,:),[col(j),'-']); hold on
   %   end
   %   title(sprintf('Display %i',d));
   %   set(gca,'TickDir','out');
   %   xlabel('Wavelength (nm)');
   %   ylabel('Power');
   %   box off
   %end
   %legend({'Primary 1' 'Primary 2' 'Primary 3'});

   %% Display matching chromatic spatial patterns
   figure('Name','Color matching');
   plotTitles = {'Display 1' 'Photopigment' 'Display 2'};
   % only plotting middle 80% of curves to avoid assumed circular edges
   idx = ceil(p.dx*0.1):floor(p.dx*0.9);
   x = p.space(idx);

   for c = 1:size(p.cones,1)
      subplot(3,3,(c-1)*3+1);
      plot(x,p.im(c,idx),[col(c),'-']);
      set(gca,'TickDir','out');
      if c==1, title('Display 1'), end
      if c==3, xlabel('Position (deg)'), end
      set(gca,'XLim',[min(p.space) max(p.space)]);
      set(gca,'YLim',[-0.05 1.05]);
      box off

      % Photopigment
      subplot(3,3,(c-1)*3+2);
      plot(x,pigSpace(c,idx),[col(c),'-']);
      set(gca,'TickDir','out');
      if c==1, title('Photopigment'), end
      if c==3, xlabel('Position (deg)'), end
      set(gca,'XLim',[min(p.space) max(p.space)]);
      box off

      % Display 2
      subplot(3,3,(c-1)*3+3);
      plot(x,pigSpace2(c,idx),[col(c),'-']);
      set(gca,'TickDir','out');
      if c==1, title('Display 2'), end
      if c==3, xlabel('Position (deg)'), end
      set(gca,'XLim',[min(p.space) max(p.space)]);
      box off
   end

end
