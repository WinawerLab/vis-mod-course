function [high_freq, low_freq,blurImg_hf, blurImg_lf] = img_blur_cones()
%% blur an image given what we know about photoreceptor density
% account for differences in density that exist in different retinal locations
%% create gratings & mark diff retinal locations in red
%clear all; clc; close all;
% top quadrant - superior
% bottom - inferior
% left - nasal
% right - temporal 
% make vertical grating of sf of 8 & 1 and 11 degs right and left of center
%remove the last row/column to make it symmetrical
[high_freq] = mk_grating(8,90,12); high_freq = high_freq(1:480,1:480);
[low_freq] = mk_grating(1,90,12); low_freq = low_freq(1:480,1:480);
%% load curcio data and normalize
dataLoad_all;
% normalize & hand pick first 12 degrees (from 10~12deg due to small
% differences in x-values for the datasets)
eccentricity = xt_deg_c(1:15);
temporal = yt_deg_c./norm(yt_deg_c); temporal = temporal(1:15);
nasal = yn_deg_c./norm(yn_deg_c); nasal = nasal(1:15);
superior = ys_deg_c./norm(ys_deg_c); superior = superior(1:15);
inferior = yi_deg_c./norm(yi_deg_c); inferior = inferior(1:15);
figure; 
loglog(eccentricity,temporal,'r',eccentricity,nasal,'g',eccentricity,superior,'b',eccentricity,inferior,'k');
%clear useless stuff
clearvars -except eccentricity temporal nasal superior inferior high_freq low_freq
%% create filters given what we know about cones
%define filters
temporal_filt = repmat(temporal,[480 32]); 
nasal_filt = repmat(nasal,[480 32]);
inferior_filt = repmat(inferior,[480 32]);
superior_filt = repmat(superior,[480 32]);
%% create a mask for each retinal location (segment into 4 parts)
close all;
rows = 480;columns = 480;
%superior 
xCoords = [0 240 480]; yCoords = [0 240 0];
mask = poly2mask(xCoords, yCoords, rows, columns);
%only take what is a one (mask is superior retina area)
superior_ret = superior_filt.*mask;
%inferior
xCoords = [0 240 480]; yCoords = [480 240 480];
mask = poly2mask(xCoords, yCoords, rows, columns);
inferior_ret = inferior_filt.*mask;
%nasal
xCoords = [0 240 0]; yCoords = [0 240 480];
mask = poly2mask(xCoords, yCoords, rows, columns);
nasal_ret = nasal_filt.*mask;
%temporal
xCoords = [480 240 480]; yCoords = [0 240 480];
mask = poly2mask(xCoords, yCoords, rows, columns);
temporal_ret = temporal_filt.*mask;
%show all separate filters
figure;
subplot(2,2,1);imshow(superior_ret); colormap('gray'); axis square; title('superior');
subplot(2,2,2);imshow(inferior_ret); colormap('gray'); axis square; title('inferior');
subplot(2,2,3);imshow(nasal_ret); axis square; title('nasal');
subplot(2,2,4);imshow(temporal_ret); colormap('gray'); axis square; title('temporal');
%put filters together
figure;
loc_filt = superior_ret+inferior_ret+nasal_ret+temporal_ret;
imagesc(loc_filt); colormap('gray'); colorbar; axis square; 
title('Full filter across all retinal locations');
%% Apply location filter to gratings
%filter index
figure; colormap('gray');
blurImg_hf = conv2(loc_filt,high_freq,'same');
blurImg_lf = conv2(loc_filt,low_freq,'same');
%% plot all
subplot(2,2,1); imagesc(high_freq); colormap('gray'); axis square; 
set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480]);
set(gca,'xticklabel',[-12,-6, 0, 6, 12]); xlabel('deg');
title('high freq grating');
subplot(2,2,2); imagesc(low_freq); axis square;
set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480]);
set(gca,'xticklabel',[-12,-6, 0, 6, 12]); xlabel('deg');
title('low freq grating');
subplot(2,2,3);imagesc(blurImg_hf); axis square; 
set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480]);
set(gca,'xticklabel',[-12,-6, 0, 6, 12]); xlabel('deg');
% set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480 600 720 840 959]);
% set(gca,'xticklabel',[-24,-18, -12, -8, 0, 8, 12, 18, 24]); xlabel('deg');
title('high sf after blur');
subplot(2,2,4);imagesc(blurImg_lf); axis square; 
set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480]);
set(gca,'xticklabel',[-12,-6, 0, 6, 12]); xlabel('deg');
% set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480 600 720 840 959]);
% set(gca,'xticklabel',[-24,-18, -12, -8, 0, 8, 12, 18, 24]); xlabel('deg');
title('low sf after blur');
%so that I can it in other functions
close all;

