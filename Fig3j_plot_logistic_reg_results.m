cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%% Parameters

Regions = {'ACAd','ILA','ORBm','PL'};
Genotypes = {'WT','VGlut2','NPY','Esr'};
load('/Users/pierre/Documents/MATLAB/neuropixelPFC/Matlab/Utilities/Colormaps/cmap_genotype2.mat')

%% load
 m_b12 = readNPY('/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/logistic_regression/decoding_acc_bl1-2_mean.npy');
 sd_b12 = readNPY('/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/logistic_regression/decoding_acc_bl1-2_sd.npy');
 
  m_b13 = readNPY('/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/logistic_regression/decoding_acc_bl1-3_mean.npy');
 sd_b13 = readNPY('/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/logistic_regression/decoding_acc_bl1-3_sd.npy');
 
  m_b14 = readNPY('/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/logistic_regression/decoding_acc_bl1-4_mean.npy');
 sd_b14 = readNPY('/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/logistic_regression/decoding_acc_bl1-4_sd.npy');
 
 m_b12(m_b12==0) = NaN; sd_b12(sd_b12==0) = NaN; 
  m_b13(m_b13==0) = NaN; sd_b13(sd_b12==0) = NaN;
   m_b14(m_b14==0) = NaN; sd_b14(sd_b12==0) = NaN;
   
 %% plot
 % block 1 vs 2
 f1 = figure;
 cnt = 0;
 for R = 1 : numel(Regions)
 subplot(1,4,1+cnt)
 errorbar(m_b12(1+cnt,:),sd_b12(1+cnt,:),'-o','Color',cmap_genotype(1,:))
 hold on
  errorbar(m_b12(5+cnt,:),sd_b12(5+cnt,:),'-o','Color',cmap_genotype(4,:))
   errorbar(m_b12(9+cnt,:),sd_b12(9+cnt,:),'-o','Color',cmap_genotype(3,:))
   errorbar(m_b12(13+cnt,:),sd_b12(13+cnt,:),'-o','Color',cmap_genotype(2,:)) 
 ylim([0.4 1])
 ylabel('accuracy')
 line([0 30],[0.5 0.5],'Color','k','LineStyle','--')
        title(Regions{1+cnt})
        set(gca,'XTick',[0 10 20 30])
        xticklabels([0 50 100 150])
        xlabel('number of units')
        cnt = cnt + 1;
 end
 
 sgtitle('Block 1 vs 2')
 set(gcf,'units','points','position',[0,800,1600,500])
 
 % block 1 vs 3
  f2 = figure;
  cnt = 0;
 for R = 1 : numel(Regions)
 subplot(1,4,1+cnt)
 errorbar(m_b13(1+cnt,:),sd_b13(1+cnt,:),'-o','Color',cmap_genotype(1,:))
 hold on
  errorbar(m_b13(5+cnt,:),sd_b13(5+cnt,:),'-o','Color',cmap_genotype(4,:))
   errorbar(m_b13(9+cnt,:),sd_b13(9+cnt,:),'-o','Color',cmap_genotype(3,:))
   errorbar(m_b13(13+cnt,:),sd_b13(13+cnt,:),'-o','Color',cmap_genotype(2,:)) 
 ylim([0.4 1])
 ylabel('accuracy')
 line([0 30],[0.5 0.5],'Color','k','LineStyle','--')
        title(Regions{1+cnt})
        set(gca,'XTick',[0 10 20 30])
        xticklabels([0 50 100 150])
        xlabel('number of units')
        cnt = cnt + 1;
 end
 
 sgtitle('Block 1 vs 3')
 set(gcf,'units','points','position',[0,800,1600,500])
 
 % block 1 vs 4
  f3 = figure;
  cnt = 0;
 for R = 1 : numel(Regions)
 subplot(1,4,1+cnt)
 errorbar(m_b14(1+cnt,:),sd_b14(1+cnt,:),'-o','Color',cmap_genotype(1,:))
 hold on
  errorbar(m_b14(5+cnt,:),sd_b14(5+cnt,:),'-o','Color',cmap_genotype(4,:))
   errorbar(m_b14(9+cnt,:),sd_b14(9+cnt,:),'-o','Color',cmap_genotype(3,:))
   errorbar(m_b14(13+cnt,:),sd_b14(13+cnt,:),'-o','Color',cmap_genotype(2,:)) 
 ylim([0.4 1])
 ylabel('accuracy')
 line([0 30],[0.5 0.5],'Color','k','LineStyle','--')
        title(Regions{1+cnt})
        set(gca,'XTick',[0 10 20 30])
        xticklabels([0 50 100 150])
        xlabel('number of units')
        cnt = cnt + 1;
 end
 
 sgtitle('Block 1 vs 4')
 set(gcf,'units','points','position',[0,800,1600,500])
 
 %% colormaps
 [cmap] = cbrewer('seq','Greens',64,'spline');
 cmap(57:64,:) = [];
  % block 1 vs 2
 f4 = figure;
 cnt = 0;
 for R = 1 : numel(Regions)
 subplot(1,4,1+cnt)
 A = nanmean(m_b12(5+cnt,:)-m_b12(1+cnt,:),2);
 B = nanmean(m_b12(9+cnt,:)-m_b12(1+cnt,:),2);
 C = nanmean(m_b12(13+cnt,:)-m_b12(1+cnt,:),2);
 
 imagesc([A B C])
 colormap(cmap)
caxis([0, 0.07])        
colorbar

 title(Regions{1+cnt})
        set(gca,'XTick',1:3)
        xticklabels({'Esr', 'NPY', 'VGlut2'})
        xlabel('Genotypes')
        cnt = cnt + 1;
 end
 
 sgtitle('Block 1 vs 2')
 set(gcf,'units','points','position',[0,800,1600,500])
 
   % block 1 vs 2
 f5 = figure;
 cnt = 0;
 for R = 1 : numel(Regions)
 subplot(1,4,1+cnt)
 A = nanmean(m_b13(5+cnt,:)-m_b13(1+cnt,:),2);
 B = nanmean(m_b13(9+cnt,:)-m_b13(1+cnt,:),2);
 C = nanmean(m_b13(13+cnt,:)-m_b13(1+cnt,:),2);
 
 imagesc([A B C])
 colormap(cmap)
caxis([0, 0.07])        
colorbar

 title(Regions{1+cnt})
        set(gca,'XTick',1:3)
        xticklabels({'Esr', 'NPY', 'VGlut2'})
        xlabel('Genotypes')
        cnt = cnt + 1;
 end
 
 sgtitle('Block 1 vs 3')
 set(gcf,'units','points','position',[0,800,1600,500])
 
   % block 1 vs 2
 f6 = figure;
 cnt = 0;
 for R = 1 : numel(Regions)
 subplot(1,4,1+cnt)
 A = nanmean(m_b14(5+cnt,:)-m_b14(1+cnt,:),2);
 B = nanmean(m_b14(9+cnt,:)-m_b14(1+cnt,:),2);
 C = nanmean(m_b14(13+cnt,:)-m_b14(1+cnt,:),2);
 
 imagesc([A B C])
 colormap(cmap)
caxis([0, 0.07])        
colorbar

 title(Regions{1+cnt})
        set(gca,'XTick',1:3)
        xticklabels({'Esr', 'NPY', 'VGlut2'})
        xlabel('Genotypes')
        cnt = cnt + 1;
 end
 
 sgtitle('Block 1 vs 4')
 set(gcf,'units','points','position',[0,800,1600,500])
        