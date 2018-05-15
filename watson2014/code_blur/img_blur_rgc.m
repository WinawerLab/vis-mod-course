function [blurImg_hf_rgcd, blurImg_lf_rgcd] = img_blur_rgc()
%% blur an image given what we know about density of rgc receptive fields
% account for differences in density that exist in different retinal locations
clear all; close all; clc;
%% re-create watson (2014) math to get rgc rf density  & load cone blur + gratings
[ogecc, fittedVals] = watson14;
[high_freq, low_freq,blurImg_hf, blurImg_lf] = img_blur_cones;
% take first ~11 deg
eccentricity = ogecc{1}; eccentricity = eccentricity(1:15);
% take first ~11 deg for all retinal locations
temporal = fittedVals{1}./norm(fittedVals{1}); temporal = temporal(1:15);
nasal = fittedVals{2}./norm(fittedVals{2}); nasal = nasal(1:15);
superior = fittedVals{3}./norm(fittedVals{3}); superior = superior(1:15);
inferior = fittedVals{4}./norm(fittedVals{4}); inferior = inferior(1:15);
figure; 
loglog(eccentricity,temporal,'r',eccentricity,nasal,'g',eccentricity,superior,'b',eccentricity,inferior,'k');
title('watson fits used to create the 4 filters');
%% create filters given what we know about retinal ganglion cell densities
%define filters
temporal_filt = repmat(temporal',[480 32]); 
nasal_filt = repmat(nasal',[480 32]);
inferior_filt = repmat(inferior',[480 32]);
superior_filt = repmat(superior',[480 32]);
%% create a mask for each retinal location (segment into 4 parts)
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
subplot(2,2,3);imshow(nasal_ret); colormap('gray'); axis square; title('nasal');
subplot(2,2,4);imshow(temporal_ret); colormap('gray'); axis square; title('temporal');
%put filters together
figure;
loc_filt = superior_ret+inferior_ret+nasal_ret+temporal_ret;
imagesc(loc_filt); colormap('gray'); axis square; 
title('Full filter across all retinal locations');
%% apply new rgc rf filters
figure; colormap('gray');
blurImg_hf_rgcd = conv2(loc_filt,blurImg_hf,'same');
blurImg_lf_rgcd = conv2(loc_filt,blurImg_lf,'same');
%% plot all
subplot(2,2,1); imagesc(high_freq); colormap('gray'); axis square; 
set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480]);
set(gca,'xticklabel',[-12,-6, 0, 6, 12]); xlabel('deg');
title('high freq grating');
subplot(2,2,2); imagesc(low_freq); axis square;
set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480]);
set(gca,'xticklabel',[-12,-6, 0, 6, 12]); xlabel('deg');
title('low freq grating');
subplot(2,2,3);imagesc(blurImg_hf_rgcd); axis square; 
set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480]);
set(gca,'xticklabel',[-12,-6, 0, 6, 12]); xlabel('deg');
% set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480 600 720 840 959]);
% set(gca,'xticklabel',[-24,-18, -12, -8, 0, 8, 12, 18, 24]); xlabel('deg');
title('high sf after blur');
subplot(2,2,4);imagesc(blurImg_lf_rgcd); axis square; 
set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480]);
set(gca,'xticklabel',[-12,-6, 0, 6, 12]); xlabel('deg');
% set(gca,'yticklabel',[]);  set(gca,'xtick',[1 120 240 360 480 600 720 840 959]);
% set(gca,'xticklabel',[-24,-18, -12, -8, 0, 8, 12, 18, 24]); xlabel('deg');
title('low sf after blur');


