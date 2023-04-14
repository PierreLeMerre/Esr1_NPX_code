%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that plots the activity modes per region from processed data.
% Inspired and adapted from Allen et al., Science, 2019
%
% Written by pielem Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%%
Regions = {'ACAd','PL','ILA','ORBm'};
Genotypes = {'WT','VGlut2','NPY','Esr'};

% activity bin parameters
bin_sz = 0.1; % in seconds
pre = -2; % in seconds
post = 5; % in seconds
time_window = post-pre; % in seconds

% preallocate
load('/Users/pielem/Documents/MATLAB/Esr1_NPX_code/analysis/activity_modes/PL_activity_modes5.mat')
sb2 = nan(16,size(state.b2,2));
cc = nan(16,size(state.b2,2));
av = nan(16,size(state.b2,2));

t = linspace(pre,post,round(time_window/bin_sz));
cnt = 1;
for g = 1:numel(Genotypes)
    
    for i = 1:numel(Regions)
        region = Regions{i};
        load(['/Users/pielem/Documents/MATLAB/Esr1_NPX_code/analysis/activity_modes/' region '_activity_modes5.mat'])
        sb2(cnt,:) = state.b2(g,:)-state.b1(g,:);
        cc(cnt,:) = cue.b2(g,:)-cue.b1(g,:);
        av(cnt,:) = avs.b2(g,:)-avs.b1(g,:);
        cnt = cnt + 1;
    end
end

figure;
cnt = 0;
cnt2 = 0;
for g = 1:numel(Genotypes)
    
    subplot(4,3,1+cnt)
    imagesc('XData',t,'CData',flip(cc(1+cnt2:4+cnt2,:)))    
    caxis([-0.15, 0.15])
    xlim([-1 2])
    set(gca,'YTick',[1 2 3 4]);
    set(gca,'YTickLabel',flip(Regions));
    title('Cue')
    
    subplot(4,3,2+cnt)
    imagesc('XData',t,'CData',flip(av(1+cnt2:4+cnt2,:)))
    xlim([-1 2])
    caxis([-0.15, 0.15])
    title('Aversive signal')
    
    subplot(4,3,3+cnt)
    imagesc('XData',t,'CData',flip(sb2(1+cnt2:4+cnt2,:)))
    xlim([-1 2])
    caxis([-0.15, 0.15])
    title('State')
    cnt = cnt + 3;
    cnt2 = cnt2 + 4;
    
    colormap(ametrine)

    
end

sgtitle('Activity mode: diff Block 2 - block 1')
set(gcf,'units','points','position',[500,370,1500,550])

%% State mode color map
%Esr
figure;
load('/Users/pielem/OneDrive - Karolinska Institutet/Mac/Documents/MATLAB/Esr1_NPX_code/utilities/Colormaps/purple.mat')
state_mean = nan(1,4);
for i = 1:4   
    state_mean(i) = mean(sb2(12+i,:),2);
end    
imagesc(state_mean)
colormap(purple_cmap)
caxis([0, 0.15])        
colorbar
title('Esr')

%WT
figure;
state_mean = nan(1,4);
for i = 1:4   
    state_mean(i) = mean(sb2(i,:),2);
end    
imagesc(state_mean)
colormap(purple_cmap)
caxis([0, 0.15])        
colorbar
title('WT')

%VGlut2
figure;
state_mean = nan(1,4);
for i = 1:4   
    state_mean(i) = mean(sb2(4+i,:),2);
end    
imagesc(state_mean)
colormap(purple_cmap)
caxis([0, 0.15])        
colorbar
title('VGlut2')

%NPY
figure;
state_mean = nan(1,4);
for i = 1:4   
    state_mean(i) = mean(sb2(8+i,:),2);
end    
imagesc(state_mean)
colormap(purple_cmap)
caxis([0, 0.15])        
colorbar
title('NPY')
