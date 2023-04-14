%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the location of the modulated units in the ARA.
%
% Written by Pierre Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%% Load

Tasks = 'Aversion'; %,'Detection','Opto'

genotype2plot = 'WT';

% Your path to the Database
Path2Data = '/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/Aversion/';
D = dir([Path2Data '*.nwb']);

% Allen Brain Atlas regions to count the units
Regions = {'ACAd','PL','ILA','ORBm'};

% Your path to analysis folder
BasePath = '/Users/pielem/Documents/MATLAB/';
Path2Ana = [BasePath 'Esr1_NPX_code/analysis/'];


MODUL_s = [];
REG_s = [];

MODUL_o = [];
REG_o = [];

MODUL_s2 = [];
REG_s2 = [];

MODUL_a = [];
REG_a = [];

ISREGION = [];

   

            %Chose genotype
            Genotype = genotype2plot;
            if strcmp(Genotype,'NPY')
                f_start = 1;
                f_stop = 5;
            elseif strcmp(Genotype,'VGlut2')
                f_start = 6;
                f_stop = 10;
            elseif strcmp(Genotype,'WT')
                f_start = 11;
                f_stop = 15;
            elseif strcmp(Genotype,'Esr-retro')
                f_start = 16;
                f_stop = 20;
            elseif strcmp(Genotype,'Esr-antero')
                f_start = 21;
                f_stop = 25;
            elseif strcmp(Genotype,'Esr')
                f_start = 16;
                f_stop = 25;
            elseif strcmp(Genotype,'All')
                f_start = 1;
                f_stop = 25;
            end


    CNT = 1;
    M_modul_s = nan(1,f_stop-f_start+1);
    M_modul_o = nan(1,f_stop-f_start+1);
    M_modul_s2 = nan(1,f_stop-f_start+1);
    M_modul_a = nan(1,f_stop-f_start+1);
    
for f = f_start : f_stop % Loop through mice
    
    % Read nwb file
    nwb = nwbRead([Path2Data D(f).name]);
    % get id and recording duration from meta file
    id = strsplit(nwb.general_session_id,'_');
    
    % Get PFC regions
    ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();
    
    % Get unit main channel
    main_ch = nwb.units.electrodes.data.load();
    
    isregion = zeros(numel(Regions),numel(main_ch));
    % loop over regions
    for R = 1:numel(Regions)
        region=char(Regions{R});
        % loop over units
        for i = 1 : numel(main_ch)
            % Make logical region vector
            if (strfind(ede_regions{main_ch(i)+1,:},region)) > 0 % main_ch is 0 indexed!
                isregion(R,i) = 1;
            end
        end
    end
    ISREGION = [ISREGION isregion];
    
    disp(['Get GLM-fitted units for mouse ' id{1} ', session ' id{2}])
    
    % load GLM units
    load([Path2Ana 'GLM_fitted/' id{1} '_' id{2} '_tmu.mat'])

    % load GLM tuning scores and p_value boleans
    load([Path2Ana 'GLM_fitted/' id{1} '_' id{2} '_tuning_scores.mat'])
    load([Path2Ana 'GLM_fitted/' id{1} '_' id{2} '_tuned_bol.mat'])
    
    
    % with all fitted units for each
    MODUL_s = [MODUL_s; tun_scores(1,:)'];
    temp_s = ede_regions(main_ch(tmu+1));
    REG_s = [REG_s; temp_s(tun_bol(1,:)==1 | tun_bol(1,:)==-1)];
    
    m1 = tun_scores(1,tun_scores(1,:)>=0);
    m2 = abs(tun_scores(1,tun_scores(1,:)<0));
    M_modul_s(CNT) = nanmean([m1 m2]);
    
    % with all fitted units for each
    MODUL_o = [MODUL_o; tun_scores(2,:)'];
    temp_o = ede_regions(main_ch(tmu+1));
    REG_o = [REG_o; temp_s(tun_bol(2,:)==1 | tun_bol(2,:)==-1)];
    m1 = tun_scores(2,tun_scores(2,:)>=0);
    %m1 = [];
    m2 = abs(tun_scores(2,tun_scores(2,:)<0));
    M_modul_o(CNT) = nanmean([m1 m2]);
    
    % with all fitted units for each
    MODUL_s2 = [MODUL_s2; tun_scores(3,:)'];
    temp_s2 = ede_regions(main_ch(tmu+1));
    REG_s2 = [REG_s2; temp_s(tun_bol(3,:)==1 | tun_bol(3,:)==-1)];
    m1 = tun_scores(3,tun_scores(3,:)>=0);
    m2 = abs(tun_scores(3,tun_scores(3,:)<0));
    M_modul_s2(CNT) = nanmean([m1 m2]);
    
    % with all fitted units for each
    MODUL_a = [MODUL_a; tun_scores(4,:)'];
    temp_a = ede_regions(main_ch(tmu+1));
    REG_a = [REG_a; temp_s(tun_bol(4,:)==1 | tun_bol(4,:)==-1)];
    m1 = tun_scores(4,tun_scores(4,:)>=0);
    m2 = abs(tun_scores(4,tun_scores(4,:)<0));
    m2(m2>30) = []; % one air puff outlier to remove
    M_modul_a(CNT) = nanmean([m1 m2]);
    
    CNT = CNT + 1;
    
end % end of mouse loop



%% PLot Positive

[cmap] = cbrewer('seq', 'Reds',64,'spline');
% sound
f1 = figure;
subplot(1,4,1)
mod_idx_s_all = nan(numel(Regions),numel(REG_s));
% Loop through regions
for R = 1:numel(Regions)
    region=char(Regions{R});
    for i = 1 : numel(REG_s)
        r = char(REG_s{i});
        %all units
        if (strfind(r,region))>0 & MODUL_s(i)>0 % main_ch is 0 indexed!
            mod_idx_s_all(R,i) = MODUL_s(i);
        end
    end
end % region loop
mod_idx_s_all = mod_idx_s_all';
imagesc(nanmean(mod_idx_s_all))
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'mean tuning score';
caxis([0, 8])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Sound')


% opto
subplot(1,4,2)
mod_idx_o_all = nan(numel(Regions),numel(REG_o));
% Loop through regions
for R = 1:numel(Regions)
    region=char(Regions{R});
    for i = 1 : numel(REG_o)
        r = char(REG_o{i});
        if (strfind(r,region)) > 0 & MODUL_o(i)>0 % main_ch is 0 indexed!
            mod_idx_o_all(R,i) = MODUL_o(i);
        end
    end
end % region loop
mod_idx_o_all = mod_idx_o_all';
imagesc(nanmean(mod_idx_o_all))
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'mean tuning score';
caxis([0, 8])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Opto')


% sound2
subplot(1,4,3)
mod_idx_s2_all = nan(numel(Regions),numel(REG_s2));
% Loop through regions
for R = 1:numel(Regions)
    region=char(Regions{R});
    for i = 1 : numel(REG_s2)
        r = char(REG_s2{i});
        if (strfind(r,region))>0 & MODUL_s2(i)>0 % main_ch is 0 indexed!
            mod_idx_s2_all(R,i) = MODUL_s2(i);
        end
    end
end % region loop
mod_idx_s2_all = mod_idx_s2_all';
imagesc(nanmean(mod_idx_s2_all))
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'mean tuning score';
caxis([0, 8])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Sound')


% air puff
subplot(1,4,4)
mod_idx_a_all = nan(numel(Regions),numel(REG_a));
% Loop through regions
for R = 1:numel(Regions)
    region=char(Regions{R});
    for i = 1 : numel(REG_a)
        r = char(REG_a{i});
        if (strfind(r,region)) > 0 & MODUL_a(i)>0% main_ch is 0 indexed!
            mod_idx_a_all(R,i) = MODUL_a(i);
        end
    end
end % region loop
mod_idx_a_all = mod_idx_a_all';
imagesc(nanmean(mod_idx_a_all))
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'mean tuning score';
caxis([0, 8])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Air puff')

sgtitle(Genotype)
set(gcf,'units','points','position',[0,2000,1500,250])

% PLot negative
[cmap] = cbrewer('seq', 'Blues',64,'spline');
% sound
f2 = figure;
subplot(1,4,1)
mod_idx_s_neg_all = nan(numel(Regions),numel(REG_s));
% Loop through regions
for R = 1:numel(Regions)
    region=char(Regions{R});
    for i = 1 : numel(REG_s)
        r = char(REG_s{i});
        %all units
        if (strfind(r,region))>0 & MODUL_s(i)<0 % main_ch is 0 indexed!
            mod_idx_s_neg_all(R,i) = abs(MODUL_s(i));
        end
    end
end % region loop
mod_idx_s_neg_all = mod_idx_s_neg_all';
imagesc(nanmean(mod_idx_s_neg_all))
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'mean tuning score';
caxis([0, 8])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Sound')

% opto
subplot(1,4,2)
mod_idx_o_neg_all = nan(numel(Regions),numel(REG_o));
% Loop through regions
for R = 1:numel(Regions)
    region=char(Regions{R});
    for i = 1 : numel(REG_o)
        r = char(REG_o{i});
        if (strfind(r,region)) > 0 & MODUL_o(i)<0 % main_ch is 0 indexed!
            mod_idx_o_neg_all(R,i) = abs(MODUL_o(i));
        end
    end
end % region loop
mod_idx_o_neg_all = mod_idx_o_neg_all';
imagesc(nanmean(mod_idx_o_neg_all))
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'mean tuning score';
caxis([0, 8])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Opto')

% sound
subplot(1,4,3)
mod_idx_s2_neg_all = nan(numel(Regions),numel(REG_s2));
% Loop through regions
for R = 1:numel(Regions)
    region=char(Regions{R});
    for i = 1 : numel(REG_s2)
        r = char(REG_s2{i});
        if (strfind(r,region))>0 & MODUL_s2(i)<0 % main_ch is 0 indexed!
            mod_idx_s2_neg_all(R,i) = abs(MODUL_s2(i));
        end
    end
end % region loop
mod_idx_s2_neg_all = mod_idx_s2_neg_all';
imagesc(nanmean(mod_idx_s2_neg_all))
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'mean tuning score';
caxis([0, 8])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Sound')

% air puff
subplot(1,4,4)
mod_idx_a_neg_all = nan(numel(Regions),numel(REG_a));
% Loop through regions
for R = 1:numel(Regions)
    region=char(Regions{R});
    for i = 1 : numel(REG_a)
        r = char(REG_a{i});
        if (strfind(r,region)) > 0 & MODUL_a(i)<0 & MODUL_a(i)>-30% main_ch is 0 indexed! %Remove one outilier here
            mod_idx_a_neg_all(R,i) = abs(MODUL_a(i));
        end
    end
end % region loop
mod_idx_a_neg_all = mod_idx_a_neg_all';
imagesc(nanmean(mod_idx_a_neg_all))
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'mean tuning score';
caxis([0, 8])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Air puff')


sgtitle(Genotype)
set(gcf,'units','points','position',[0,775,1500,250])


 %% count units
[cmap] = cbrewer('seq', 'Reds',64,'spline');
f3 = figure;
subplot(1,4,1)
sum_s = sum(~isnan(mod_idx_s_all)) + sum(~isnan(mod_idx_s_neg_all));
sum_s_all = sum(ISREGION,2)';
imagesc(sum_s./sum_s_all);
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'fraction of tuned units';
caxis([0, .3])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Sound')

subplot(1,4,2)
sum_o = sum(~isnan(mod_idx_o_all)) + sum(~isnan(mod_idx_o_neg_all));
sum_o_all = sum(ISREGION,2)';
imagesc(sum_o./sum_o_all);
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'fraction of tuned units';
caxis([0, .25])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Opto')

subplot(1,4,3)
sum_s2 = sum(~isnan(mod_idx_s2_all)) + sum(~isnan(mod_idx_s2_neg_all));
sum_s2_all = sum(ISREGION,2)';
imagesc(sum_s2./sum_s2_all);
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'fraction of tuned units';
caxis([0, .25])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Sound2')

subplot(1,4,4)
sum_a = sum(~isnan(mod_idx_a_all)) + sum(~isnan(mod_idx_a_neg_all));
sum_a_all = sum(ISREGION,2)';
imagesc(sum_a./sum_a_all);
colormap(cmap)
c=colorbar('southoutside');
c.Label.String = 'fraction of tuned units';
caxis([0, .25])
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',Regions);
set(gca,'YTick',[]);
title('Airpuff')

sgtitle(Genotype)
set(gcf,'units','points','position',[0,475,1500,250])

