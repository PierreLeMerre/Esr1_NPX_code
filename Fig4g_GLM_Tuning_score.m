%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find beta weigths of the GLM significanlty different from randomized
% sessions
%
%
% Written by pierre Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%% Load

% Task
Task = 'Aversion';

% Your path to the Database
Path2Data = '/Volumes/T7/NPX_Database/mPFC/Aversion/';
D = dir([Path2Data '*.nwb']);

% Your path to analysis folder
BasePath = '/Users/pierre/Documents/MATLAB/';
Path2Ana = [BasePath 'Esr1_NPX_code/analysis/'];

% to save csv 
dest = '/Users/pierre/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Pierre_Shared/LHA-LHb-PFC/Nature_Neuroscience_Revisions/tuning_scores_csv/';

% Load colormaps
load([BasePath 'Esr1_NPX_code/utilities/Colormaps/sd_colormap.mat'])


% Randomization
nboot = 1000; % nb of surrogate data

% Allen Brain Atlas regions to count the units
Regions = {'ACAd','PL','ILA','ORBm'};

% predictors
predictors = {'Sound','Opto','Sound2','Air puff'};
%predictors = {'Sound'};
%Chose genotype
Genotype = 'Esr';
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
    f_start = 22;
    f_stop = 26;
elseif strcmp(Genotype,'Esr-antero')
    f_start = 27;
    f_stop = 31;
elseif strcmp(Genotype,'Esr')
    f_start = 22;
    f_stop = 31;
elseif strcmp(Genotype,'All')
    f_start = 1;
    f_stop = 31;
end


T = [];
S = [];
TS = [];
RS = [];
FS = [];
CNT = 0;

for f = f_start : f_stop
    
    if strcmp(Genotype,'All') && 15<f && f<22
        
    else
        
        % Read nwb file
        nwb = nwbRead([Path2Data D(f).name]);
        % get id and recording duration from meta file
        id = strsplit(nwb.general_session_id,'_');
        id2 = strsplit(id{2},'-');
        disp(['Reading ' id])
        % Get PFC regions
        ede_regions = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();
        % Get unit main channel
        main_ch = nwb.units.electrodes.data.load();
        isregion_mtrx = zeros(numel(main_ch),numel(Regions));
        
        %FS and RS
        fs_bol = nwb.units.vectordata.get('FS').data.load;
        rs_bol = nwb.units.vectordata.get('RS').data.load;
        
        % loop over regions
        for R = 1:numel(Regions)
            region=char(Regions{R});
            isregion = zeros(1,numel(main_ch));
            
            % loop over units
            for i = 1 : numel(main_ch)
                
                % Make logical region vector
                if (strfind(char(ede_regions(main_ch(i)+1,:)),region)) > 0 % main_ch is 0 indexed!
                    isregion(i) = 1;
                end
                
            end
            isregion_mtrx(:,R) = isregion;
        end
        
        
        load(['/Users/pielem/Documents/MATLAB/neuropixelPFC/Matlab/analysis/GLM/' nwb.general_session_id '.mat'])
        
        %classify a unit as positively tuned to a predictor if the beta coefficient form the GLM using the
        %real data is larger than in 99.99% of the circular shifted iterations or, reversely, as negatively tuned if it
        %is less than 99.99% of the shuffled iterations.
        % tuning score is the real average minus the grand average of the circular shifted iterations and
        % divided by the standard deviation of the averages of the circular shifted iterations
        tuned = nan(numel(predictors),size(beta,2));
        shuff_tuning_scores= nan(1000,numel(predictors),size(beta,2));
        tuning_scores = nan(numel(predictors),size(beta,2));
        tmu = nan(1,size(beta,2));
        rs_glm = nan(1,size(beta,2));
        fs_glm = nan(1,size(beta,2));
        % Loop through regions
        for R = 1:numel(Regions)
            % Loop through units
            for u = 1:size(beta,2)
                if isregion_mtrx(u,R) == 1
                    CNT = CNT + 1;
                    tmu(u) = u-1; % tmu is 0 indexed !!
                    rs_glm(u) = rs_bol(u);
                    fs_glm(u) = fs_bol(u);
                    % loop through predictors
                    for p = 2:numel(predictors)+1
                        if ~isnan(beta(1,u,p))
                            
                            % upper percentile
                            up_perc = prctile(beta(2:1001,u,p),99.9);
                            % lower percentile
                            low_perc = prctile(beta(2:1001,u,p),0.01);
                            if beta(1,u,p)>=up_perc
                                tuned(p-1,u) = 1;
                            elseif beta(1,u,p)<=low_perc
                                tuned(p-1,u) = -1;
                            elseif beta(1,u,p)<up_perc && beta(1,u,p)>low_perc
                                tuned(p-1,u) = 0;
                            end
                            % tuning score
                            tuning_scores(p-1,u) = (beta(1,u,p)-nanmean(beta(2:1001,u,p)))/nanstd(beta(2:1001,u,p));
                            for i = 1 : 1000
                                shuff_tuning_scores(i,p-1,u) = (beta(i+1,u,p)-nanmean(beta(2:1001,u,p)))/nanstd(beta(2:1001,u,p));
                            end
                        end
                    end
                end
            end
        end
        
        T = cat(2,T,tuned);
        S = cat(3,S,shuff_tuning_scores);
        TS = cat(2,TS,tuning_scores);
        RS = cat(2,RS,rs_glm);
        FS = cat(2,FS,fs_glm);
        
        tun_scores = tuning_scores;
        tun_bol = tuned;
        tmu(isnan(tun_scores(1,:)))=[];
        rs_glm(isnan(tun_scores(1,:)))=[];
        fs_glm(isnan(tun_scores(1,:)))=[];
        save([Path2Ana '/GLM_fitted/' nwb.general_session_id '_tmu.mat'],'tmu')
        save([Path2Ana '/GLM_fitted/' nwb.general_session_id '_rs_glm.mat'],'rs_glm')
        save([Path2Ana '/GLM_fitted/' nwb.general_session_id '_fs_glm.mat'],'fs_glm')
        tun_scores(:,isnan(tun_scores(1,:)))=[];
        tun_bol(:,isnan(tun_bol(1,:)))=[];
        save([Path2Ana '/GLM_fitted/' nwb.general_session_id '_tuning_scores.mat'],'tun_scores')
        save([Path2Ana '/GLM_fitted/' nwb.general_session_id '_tuned_bol.mat'],'tun_bol')
        
    end
end

T(:,isnan(T(1,:)))=[];
%S2 = reshape(S(:,1,:),1,[]);
%S2(isnan(S2))=[];
FS(:,isnan(TS(1,:)))=[];
RS(:,isnan(TS(1,:)))=[];
TS(:,isnan(TS(1,:)))=[];


%% Plot example tuning scores for Sound
figure;
ee = 20;
S2 = reshape(S(:,2,:),1,[]);
S2(isnan(S2))=[];
[count1, edges1] = histcounts(S2,ee*2+1,'BinLimits',[-ee ee],'BinEdges',-ee:1:ee,'Normalization','probability');
hold on
[count2, edges2] = histcounts(TS(2,:),ee*2+1,'BinLimits',[-ee ee],'BinEdges',-ee:1:ee,'Normalization','probability');
plot(edges1(1:end-1)+1/2,count1,'k','Linewidth',2)
%plot(edges2(1:end-1)+1/2,count2,'Color',[.5 .5 .5],'Linewidth',2)
%histogram(TS(1,:),ee*2+1,'FaceColor',[.5 .5 .5],'BinLimits',[-ee ee],'BinEdges',-ee:1:ee,'Normalization','probability')
ccc = histcounts(TS(2,:),ee*2+1,'BinLimits',[-ee ee],'BinEdges',-ee:1:ee,'Normalization','probability');
b = bar(edges2(1:end-1)+1/2,ccc,'FaceColor','flat');
for i = 1:size(b.CData,1)
    if i<18
        b.CData(i,:) = [0 0 .8];
    elseif i>22
        b.CData(i,:) = [.8 0 0];
    else
        b.CData(i,:) = [.5 .5 .5];
    end
end
line([0 0],[0 0.4],'Color','k','LineStyle',':')
%histogram(TS(1,T(1,:)==1),ee*2+1,'FaceColor',[.8 0 0],'BinLimits',[-ee ee],'BinEdges',-ee:1:ee,'Normalization','probability')
%histogram(TS(1,T(1,:)==-1),ee*2+1,'FaceColor',[0 0 .8],'BinLimits',[-ee ee],'BinEdges',-ee:1:ee,'Normalization','probability')
xlabel('tuning score')
ylabel('density')


title('Sound tuning')


%% Tunings to all predictors
figure;
for p = 1:numel(predictors)
    subplot(1,numel(predictors),p)
    ee = 20;
    histogram(TS(p,:),ee*2+1,'FaceColor','k','BinLimits',[-ee ee],'BinEdges',-ee:1:ee)
    xlim([-20 20])
    hold on
    histogram(TS(p,T(p,:)==1),ee*2+1,'FaceColor','r','BinLimits',[-ee ee],'BinEdges',-ee:1:ee)
    histogram(TS(p,T(p,:)==-1),ee*2+1,'FaceColor','b','BinLimits',[-ee ee],'BinEdges',-ee:1:ee)
    xlabel('tuning score')
    ylabel('unit count')
    title(predictors{p})
end

if strcmp(Task,'Aversion')
    sgtitle(Genotype)
else
    
end


%% Pie charts of primary tunings
maxloc = zeros(1,size(TS,2));
for i = 1 : size(TS,2)
    if any(T(:,i)==1 | T(:,i)==-1)
        [~,maxloc(i)] = max(TS(:,i));
    end
end
figure;
explode = [1 1 1 1];
pie([sum(maxloc==1)/size(TS,2) sum(maxloc==2)/size(TS,2) ...
    sum(maxloc==3)/size(TS,2) sum(maxloc==4)/size(TS,2) ...
    ],explode,{'Sound 1', 'Opto', 'Sound 2', 'Airpuff'})
title(['Primary tunings ' Genotype])


%% Scatter plots

load('/Users/pierre/Documents/MATLAB/Colormaps/scatter_red_cmap.mat')
load('/Users/pierre/Documents/MATLAB/Colormaps/scatter_blue_cmap.mat')

for i = 1 : size(TS,1)
eval(['x= TS(' num2str(i) ',:);'])
eval(['c' num2str(i) ' = nan(numel(x),3);'])
for cc = 1 : numel(x)
    if (round(x(cc)/10*63)+1)>64
        eval(['c' num2str(i) '(cc,:) =  scatter_red_cmap(64,:);'])
    elseif (-1*round(x(cc)/10*63)+1)>64
        eval(['c' num2str(i) '(cc,:) =  scatter_blue_cmap(64,:);'])
    else
        if x(cc)>0
            eval(['c' num2str(i) '(cc,:) =  scatter_red_cmap(round(x(cc)/10*63)+1,:);'])
        elseif x(cc)<0
            eval(['c' num2str(i) '(cc,:) =  scatter_blue_cmap(-1*round(x(cc)/10*63)+1,:);'])
        end
    end
end
end

%figure
%imagesc(c1)
%colormap(scatter_blue_cmap)
%colorbar



f1=figure;
hold on
% scatter(TS(1,RS==1),TS(2,RS==1),[],'MarkerEdgeColor',[0 0 0.5],...
%             'MarkerFaceColor','none',...
%             'LineWidth',2,'SizeData', 200);
% scatter(TS(1,FS==1),TS(2,FS==1),[],'MarkerEdgeColor',[1 0.5 0],...
%             'MarkerFaceColor','none',...
%             'LineWidth',2,'SizeData', 200);
scatter(TS(1,:),TS(2,:),'MarkerEdgeColor','none',...
    'MarkerFaceColor',[.8 .8 .8],...
    'SizeData', 200);        
scatter(TS(1,T(1,:)==1),TS(2,T(1,:)==1),[],c1(T(1,:)==1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
scatter(TS(1,T(1,:)==-1),TS(2,T(1,:)==-1),[],c1(T(1,:)==-1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
%scatter(TS(1,T(2,:)==1),TS(2,T(2,:)==1),[],c2(T(2,:)==1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
%scatter(TS(1,T(2,:)==-1),TS(2,T(2,:)==-1),[],c2(T(2,:)==-1,:),'filled','MarkerEdgeColor','none','SizeData', 200);

%scatter(TS(1,T(1,:)==1 & T(2,:)==1),TS(2,T(1,:)==1 & T(2,:)==1),[],(c(T(1,:)==1 & T(2,:)==1,:)+c3(T(1,:)==1 & T(2,:)==1,:))/2,'filled','MarkerEdgeColor','none','SizeData', 200);

line([0 0],[-20 20],'Color','k')
line([-20 20],[0 0],'Color','k')
xlabel(['tuning score ' predictors{1}])
xlim([-20 20])
ylabel(['tuning score ' predictors{2}])
ylim([-20 20])

%sgtitle(Genotype)
set(gcf,'units','points','position',[0,2000,500,550])
%saveas(f1,['/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/GLM_Scatter_' Genotype '_sound1-opto.eps'],'epsc')


if numel(predictors)>3
% sound 2 vs airpuff
f2 = figure;

hold on
% scatter(TS(3,RS==1),TS(4,RS==1),[],'MarkerEdgeColor',[0 0 0.5],...
%             'MarkerFaceColor','none',...
%             'LineWidth',2,'SizeData', 200);
% scatter(TS(3,FS==1),TS(4,FS==1),[],'MarkerEdgeColor',[1 0.5 0],...
%             'MarkerFaceColor','none',...
%             'LineWidth',2,'SizeData', 200);
scatter(TS(3,:),TS(4,:),'MarkerEdgeColor','none',...
    'MarkerFaceColor',[.8 .8 .8],...
    'SizeData', 200);
scatter(TS(3,T(3,:)==1),TS(4,T(3,:)==1),[],c3(T(3,:)==1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
scatter(TS(3,T(3,:)==-1),TS(4,T(3,:)==-1),[],c3(T(3,:)==-1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
%scatter(TS(3,T(4,:)==1),TS(4,T(4,:)==1),[],c4(T(4,:)==1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
%scatter(TS(3,T(4,:)==-1),TS(4,T(4,:)==-1),[],c4(T(4,:)==-1,:),'filled','MarkerEdgeColor','none','SizeData', 200);

line([0 0],[-20 20],'Color','k')
line([-20 20],[0 0],'Color','k')
xlabel(['tuning score ' predictors{3}])
xlim([-10 10])
ylabel(['tuning score ' predictors{4}])
ylim([-11 11])
%sgtitle(Genotype)
set(gcf,'units','points','position',[500,2000,500,550])
%saveas(f2,['/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/GLM_Scatter_' Genotype '_sound2-airpuff.eps'],'epsc')


% opto vs Airpuff
f3 = figure;
hold on
% scatter(TS(2,RS==1),TS(4,RS==1),[],'MarkerEdgeColor',[0 0 0.5],...
%             'MarkerFaceColor','none',...
%             'LineWidth',2,'SizeData', 200);
% scatter(TS(2,FS==1),TS(4,FS==1),[],'MarkerEdgeColor',[1 0.5 0],...
%             'MarkerFaceColor','none',...
%             'LineWidth',2,'SizeData', 200);
scatter(TS(2,:),TS(4,:),'MarkerEdgeColor','none',...
    'MarkerFaceColor',[.8 .8 .8],...
    'SizeData', 200);        
scatter(TS(2,T(2,:)==1),TS(4,T(2,:)==1),[],c2(T(2,:)==1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
scatter(TS(2,T(2,:)==-1),TS(4,T(2,:)==-1),[],c2(T(2,:)==-1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
line([0 0],[-20 20],'Color','k')
line([-20 20],[0 0],'Color','k')
xlabel(['tuning score ' predictors{2}])
xlim([-20 20])
ylabel(['tuning score ' predictors{4}])
ylim([-20 20])
%sgtitle(Genotype)
set(gcf,'units','points','position',[1000,2000,500,550])

%saveas(f3,['/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/GLM_Scatter_' Genotype '_opto-airpuff.eps'],'epsc')


%Save csv data
% tuning scores
%         m = TS(2,:);
%         m2 = TS(4,:);
%         id = 1 : size(m,2);
%         sig = T(2,:);
%         values = [id' m' m2' sig'];
%         Table = array2table(values);
%         Table.Properties.VariableNames(1:4) = {'id','tuning score Opto','tuning score Air puff','significance'};
%         writetable(Table,[dest '/tuning_scores_' Genotype '.csv'],'Delimiter',',')

end

[rho, pval]= corr(TS(2,:)',TS(4,:)','type','Spearman','rows','complete');

f1=figure;
hold on
% scatter(TS(1,RS==1),TS(2,RS==1),[],'MarkerEdgeColor',[0 0 0.5],...
%             'MarkerFaceColor','none',...
%             'LineWidth',2,'SizeData', 200);
% scatter(TS(1,FS==1),TS(2,FS==1),[],'MarkerEdgeColor',[1 0.5 0],...
%             'MarkerFaceColor','none',...
%             'LineWidth',2,'SizeData', 200);
scatter(TS(1,:),TS(2,:),'MarkerEdgeColor','none',...
    'MarkerFaceColor',[.8 .8 .8],...
    'SizeData', 200);        
scatter(TS(1,T(1,:)==1),TS(2,T(1,:)==1),[],c1(T(1,:)==1,:),'filled','MarkerEdgeColor','none','SizeData', 200);
scatter(TS(1,T(1,:)==-1),TS(2,T(1,:)==-1),[],c1(T(1,:)==-1,:),'filled','MarkerEdgeColor','none','SizeData', 200);

line([0 0],[-20 20],'Color','k')
line([-20 20],[0 0],'Color','k')
xlabel(['tuning score ' predictors{1}])
xlim([-20 20])
ylabel(['tuning score ' predictors{2}])
ylim([-20 20])

%sgtitle(Genotype)
set(gcf,'units','points','position',[0,2000,500,550])
%saveas(f1,['/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/GLM_Scatter_' Genotype '_sound1-opto.eps'],'epsc')



%% Show fraction for scatter plot

figure;
b2plot = [sum(T(2,:)==1)/numel(TS(2,:)) -1*sum(T(2,:)==-1)/numel(TS(2,:))];
bar(b2plot)
ylabel('fraction of GLM fitted units')
xlabel('Opto')
ylim([-0.2 0.5])

%%
%  [cmap_blue] = cbrewer('seq','Blues',64,'spline');
%  [cmap_red] = cbrewer('seq','Reds',64,'spline');


%% Bar plots with mean tuning scores for positve and negative
% pos_table = [];
% % positive
% for i = 1 : numel(predictors)
% a = TS(i,T(i,:)==1);
% b = TS(i,T(i,:)==1 & RS==1);
% c = TS(i,T(i,:)==1 & FS==1);
% 
% p1 = ranksum(a,b)*3; p2 = ranksum(a,c)*3; p3 = ranksum(b,c)*3;
% 
% f1 = figure;
% bar([mean(a) mean(b) mean(c)],'k')
% hold on
% errorbar([mean(a) mean(b) mean(c)],[std(a) std(b) std(c)],'ko')
%  ylabel('Mean pos. tuning score')
%  ylim([0 10])
%     set(gca,'XTick',1:3)
%     set(gca,'XTickLabel',{'All', 'RS', 'FS'});
% xlabel('units')
% if i==1
% title([predictors{1} ', pvalues: ' num2str(round(p1,3)) ', ' num2str(round(p2,3)) ', ' num2str(round(p3,3))])
% %saveas(f1,[FigPath 'GLM_ts_FS_RS_' Genotype '_s'],'epsc')
% elseif i==2
% title([predictors{2} ', pvalues: ' num2str(round(p1,3)) ', ' num2str(round(p2,3)) ', ' num2str(round(p3,3))])
% %saveas(f1,[FigPath 'GLM_ts_FS_RS_' Genotype '_o'],'epsc')
% elseif i==3
% title([predictors{3} ', pvalues: ' num2str(round(p1,3)) ', ' num2str(round(p2,3)) ', ' num2str(round(p3,3))])
% %saveas(f1,[FigPath 'GLM_ts_FS_RS_' Genotype '_s2'],'epsc')
% elseif i==4
% title([predictors{4} ', pvalues: ' num2str(round(p1,3)) ', ' num2str(round(p2,3)) ', ' num2str(round(p3,3))])
% %saveas(f1,[FigPath 'GLM_ts_FS_RS_' Genotype '_a'],'epsc')
% end
% 
% pos_table(i,:) = round([mean(a) std(a) mean(b) std(b) p1 mean(c) std(c) p2 p3],3);
% 
% end
% 
% % negative
% neg_table = [];
% 
% for i = 1 : numel(predictors)
% a = TS(i,T(i,:)==-1 & TS(i,:)>=-30); %remove outlier here
% b = TS(i,T(i,:)==-1 & TS(i,:)>=-30 & RS==1);
% c = TS(i,T(i,:)==-1 & TS(i,:)>=-30 & FS==1);
% 
% 
% if ~isempty(c)
% p1 = ranksum(a,b)*3; p2 = ranksum(a,c)*3; p3 = ranksum(b,c)*3;
% else
% p1 = ranksum(a,b)*3; p2 = NaN; p3 = NaN;
% end
% 
% f2=figure;
% bar([mean(a) mean(b) mean(c)],'k')
% hold on
% errorbar([mean(a) mean(b) mean(c)],[std(a) std(b) std(c)],'ko')
%  ylabel('Mean neg. tuning score')
%   ylim([-10 0])
%     set(gca,'XTick',1:3)
%     set(gca,'XTickLabel',{'All', 'RS', 'FS'});
% xlabel('units')
% if i==1
% title([predictors{1} ', pvalues: ' num2str(round(p1,2)) ', ' num2str(round(p2,2)) ', ' num2str(round(p3,2))])
% %saveas(f2,[FigPath 'GLM_ts_FS_RS_' Genotype '_s_neg'],'epsc')
% elseif i==2
% title([predictors{2} ', pvalues: ' num2str(round(p1,2)) ', ' num2str(round(p2,2)) ', ' num2str(round(p3,2))])
% %saveas(f2,[FigPath 'GLM_ts_FS_RS_' Genotype '_o_neg'],'epsc')
% elseif i==3
% title([predictors{3} ', pvalues: ' num2str(round(p1,2)) ', ' num2str(round(p2,2)) ', ' num2str(round(p3,2))])
% %saveas(f2,[FigPath 'GLM_ts_FS_RS_' Genotype '_s2_neg'],'epsc')
% elseif i==4
% title([predictors{4} ', pvalues: ' num2str(round(p1,2)) ', ' num2str(round(p2,2)) ', ' num2str(round(p3,2))])
% %saveas(f2,[FigPath 'GLM_ts_FS_RS_' Genotype '_a_neg'],'epsc')
% end
% 
% neg_table(i,:) = round([mean(a) std(a) mean(b) std(b) p1 mean(c) std(c) p2 p3],3);
% 
% end
