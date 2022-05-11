%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the location of the modulated units in the ARA.
%
% Written by Pierre Le Merre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cls
format long
set(0, 'DefaultFigureRenderer', 'painters')

%%
Genotypes = {'WT','VGlut2','NPY','Esr'};
Stims = {'s','o','s2','a'};

for i = 1 : numel(Stims)
    Stim = Stims{i};
    array2test = nan(4,10);
    for g = 1: numel(Genotypes)
        Gen = Genotypes{g};
        load(['/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/Stats/Modul_GLMneg_' Stim '_' Gen '.mat'])
        if i == 1
            array2test(g,1:numel(M_modul_s)) = M_modul_s;
        elseif i == 2
            array2test(g,1:numel(M_modul_o)) = M_modul_o;
        elseif i == 3
            array2test(g,1:numel(M_modul_s2)) = M_modul_s2;
        elseif i == 4
            array2test(g,1:numel(M_modul_a)) = M_modul_a;
        end
    end
    array2test = array2test';
    M = nanmean(array2test);
    S = nanstd(array2test);
    figure;
    plot([.8 .8 .8 .8 .8],array2test(1:5,1),'ko')
    hold on
    plot([1.8 1.8 1.8 1.8 1.8],array2test(1:5,2),'ko')
    plot([2.8 2.8 2.8 2.8 2.8],array2test(1:5,3),'ko')
    plot([3.8 3.8 3.8 3.8 3.8 3.8 3.8 3.8 3.8 3.8],array2test(1:10,4),'ko')
    errorbar(M,S,'ko')
    xlim([.5 4.5])
    P1 = ranksum(array2test(~isnan(array2test(:,1)),1),array2test(~isnan(array2test(:,2)),2),'TAIL','right')*3;
    P2 = ranksum(array2test(~isnan(array2test(:,1)),1),array2test(~isnan(array2test(:,3)),3),'TAIL','right')*3;
    P3 = ranksum(array2test(~isnan(array2test(:,1)),1),array2test(~isnan(array2test(:,4)),4),'TAIL','right')*3;
    disp(['Stim ' Stim ])
    disp(['Mean = ' num2str(M)])
    disp(['Std = ' num2str(S)])
    disp(['P values = ' num2str(round(P1,3)), ', ' num2str(round(P2,3)) ', ' num2str(round(P3,3))])
    
end



%% Pooled units
% Genotypes = {'WT','VGlut2','NPY','Esr'};
% Stims = {'s','o','s2','a'};
% 
% for i = 1 : numel(Stims)
%     Stim = Stims{i};
%     array2test = nan(4,800);
%     for g = 1: numel(Genotypes)
%         Gen = Genotypes{g};
%         load(['/Users/pierre/OneDrive - KI.SE/Pierre_Shared/LHA-LHb-PFC/NatureSubmission/Stats/Modul_GLMpool_' Stim '_' Gen '.mat'])
%         if i == 1
%             MODUL_s(MODUL_s==0) = NaN;
%             MODUL_s = abs(MODUL_s);
%             array2test(g,1:numel(MODUL_s)) = MODUL_s;
%         elseif i == 2
%             MODUL_o(MODUL_o==0) = NaN;
%             MODUL_o = abs(MODUL_o);
%             array2test(g,1:numel(MODUL_o)) = MODUL_o;
%         elseif i == 3
%             MODUL_s2(MODUL_s2==0) = NaN;
%             MODUL_s2 = abs(MODUL_s2);
%             array2test(g,1:numel(MODUL_s2)) = MODUL_s2;
%         elseif i == 4
%             MODUL_a(MODUL_a==0) = NaN;
%             MODUL_a = abs(MODUL_a);
%             MODUL_a(MODUL_a>=30) = NaN;
%             array2test(g,1:numel(MODUL_a)) = MODUL_a;
%         end
%     end
%     
%     array2test = array2test';
%     M = nanmean(array2test);
%     S = nanstd(array2test);
%     figure;    
%     plot(.8*ones(1,800),array2test(:,1),'ko')
%     hold on
%     plot(1.8*ones(1,800),array2test(:,2),'ko')
%     plot(2.8*ones(1,800),array2test(:,3),'ko')
%     plot(3.8*ones(1,800),array2test(:,4),'ko')
%     errorbar(M,S,'ko')
%     xlim([.5 4.5])
%     P1 = ranksum(array2test(~isnan(array2test(:,1)),1),array2test(~isnan(array2test(:,2)),2))*3;
%     P2 = ranksum(array2test(~isnan(array2test(:,1)),1),array2test(~isnan(array2test(:,3)),3))*3;
%     P3 = ranksum(array2test(~isnan(array2test(:,1)),1),array2test(~isnan(array2test(:,4)),4))*3;
%     disp(['Stim ' Stim ])
%     disp(['Mean = ' num2str(M)])
%     disp(['Std = ' num2str(S)])
%     disp(['P values = ' num2str(round(P1,8)), ', ' num2str(round(P2,8)) ', ' num2str(round(P3,8))])
%     
% end