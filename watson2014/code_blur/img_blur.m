%% blur an image given what we know about photoreceptor and rgc density
% account for differences in density that exist in different retinal locations
%% create gratings & mark diff retinal locations in red
clear all; clc; close all;
% top quadrant - superior
% bottom - inferior
% left - nasal
% right - temporal 
% make vertical grating of sf of 8 & 2 and 11 degs right and left of center
[high_freq] = mk_grating(8,90,12);
[low_freq] = mk_grating(2,90,12);
%plot gratings
subplot(1,2,1); imagesc(high_freq); colormap('gray'); axis square;
line([481 0],[0 481],'color','r','linewidth',1); 
line([481 0],[481 0],'color','r','linewidth',1);
title('high freq grating');
subplot(1,2,2); imagesc(low_freq); axis square;
line([481 0],[0 481],'color','r','linewidth',1); 
line([481 0],[481 0],'color','r','linewidth',1);
title('low freq grating');
%% load curcio data and normalize
dataLoad_all;
% normalize & hand pick first 18 degrees (from 10~12deg due to small
% differences in x-values for the datasets)
temporal = yt_deg_c./norm(yt_deg_c); temporal = temporal(1:15);
nasal = yn_deg_c./norm(yn_deg_c); nasal = nasal(1:15);
superior = ys_deg_c./norm(ys_deg_c); superior = superior(1:15);
inferior = yi_deg_c./norm(yi_deg_c); inferior = inferior(1:15);
%clear useless stuff
clearvars -except eccentricity temporal nasal superior inferior high_freq low_freq
%% create filters given what we know about cones
%define filters
temporal_filt = repmat(temporal,[481 32]);
nasal_filt = repmat(nasal,[481 32]);
inferior_filt = repmat(inferior,[481 32]);
superior_filt = repmat(superior,[481 32]);
%filter index
filtNames = {'temporal_filt','nasal_filt','inferior_filt','superior_filt'};
figure; colormap('gray');
for ii = 1:4 
    blurImg_hf = conv2(eval(filtNames{ii}),high_freq,'same');
    blurImg_lf = conv2(eval(filtNames{ii}),low_freq,'same');
    subplot(1,2,1);imagesc(blurImg_hf); hold on; axis square
    line([480 0],[0 481],'color','r','linewidth',1); line([480 0],[481 0],'color','r','linewidth',1);
    subplot(1,2,2);imagesc(blurImg_lf); hold on; axis square
    line([480 0],[0 481],'color','r','linewidth',1); line([480 0],[481 0],'color','r','linewidth',1);
end

