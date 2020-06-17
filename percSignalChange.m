%% 
% This script takes the extracted ROI voxel data and shows the average 
% across groups, comparing whether there are significant differences.


clear, clc

%% Parameters

task = 'video';
group = {'P','C'};
% group = {'P'};
% group = {'C'};
smoothing = 's4';

signal_mode = 'PercentSignalChange';


%% Directories and contrasts

voxelDataSaveDir = []; % Directory where the voxel data was saved
resultDir = []; % Directory to save results of this script


if ~exist(resultDir,'dir'), mkdir(resultDir), end

contrast_names = {'SvsC','S1vsC1','S2vsC2'}; % Name of the contrasts from which the voxel data was extracted. 

%% Functional ROI comparisons

save_figure = 0;

for iCon = 1:length(contrast_names)
    
    % Loads voxel data and subject data
    % voxel data is a struct with the fields expVoxelData (voxel data for
    % the experimental condition) and conVoxelData (voxel data for the 
    % control condition). Substracting the two of them gives us the
    % contrast data.
    contrast_no = contrast_names{iCon};
    voxelDataSaveFile = dir(fullfile(voxelDataSaveDir,['ROI_data_joinROIs_' signal_mode '_' contrast_no '_*']));
    load(fullfile(voxelDataSaveFile(end).folder,voxelDataSaveFile(end).name))
    
    nCol = 5;
    nRow = ceil(length(roiLabel)/nCol);
    nFig = nCol*nRow;
            
    % Indices to subjects belonging to each group 
    grIdx = unique([S.Group]);
    cIdx = find([S.Group]==grIdx(1));
    pIdx = find([S.Group]==grIdx(2));
        
    fh = figure('PaperType','a4','PaperOrientation','landscape',...
                'PaperUnit','normalized','PaperPosition',[0 0 1 1]);
    annotation('textbox',[0.45, 0.05,0.9, 0.95],'String',contrast_names{iCon},'FitBoxToText','on')
    
    % Box plots
    for iROI = 1:length(roiLabel)
        
        fprintf('Extracting data for contrast %s - ROI:%s...\n',contrast_no,roiLabel{iROI})
        
        roiVoxelData = reshape(mean([voxelData_all(:,iROI).expVoxelData]),[],1)-reshape(mean([voxelData_all(:,iROI).conVoxelData]),[],1);
        cG = roiVoxelData(cIdx);
        pG = roiVoxelData(pIdx);
        
        [h(iCon,iROI),p(iCon,iROI),ci(iCon,iROI,:),stats(iCon,iROI)] = ttest2(pG,cG,'Alpha',0.05); % Correct Alpha due to number of comparisons              
       
        subplot(nRow,nCol,mod(iROI-1,nFig)+1)
        boxplot([pG',cG'],[zeros(1,length(pG)) ones(1,length(cG))]);
        title([roiLabel{iROI} ' - p=' num2str(p(iROI),4)],'Interpreter','none')
        set(gca,'xticklabel',{'Papa','Control'})
        ylabel('%BOLD signal change')
        
        cD(iCon,iROI,:) = cG;
        pD(iCon,iROI,:) = pG;
        
    end
    
    if save_figure
        print(fh,fullfile(resultDir,[contrast_names{iCon} '_box_' datestr(datetime('now'),'yymmddTHHMMSS') '.pdf']),'-dpdf')
        savefig(fh,fullfile(resultDir,[contrast_names{iCon} '_box_' datestr(datetime('now'),'yymmddTHHMMSS') '.fig']))
    end
    
    % Bar plots
    cD1 = squeeze(cD(iCon,:,:));
    pD1 = squeeze(pD(iCon,:,:));

    fh = figure('PaperType','a4','PaperOrientation','landscape',...
                'PaperUnit','normalized','PaperPosition',[0 0 1 1]);
%     annotation('textbox',[0.45, 0.05,0.9, 0.95],'String',contrast_names{iCon},'FitBoxToText','on')
        
    cD1 = filloutliers(cD1,'clip');
    pD1 = filloutliers(pD1,'clip');
    
    mRes = [mean(cD1,2),mean(pD1,2)];
    sRes = [std(cD1,[],2)/sqrt(size(cD,3)),std(pD1,[],2)/sqrt(size(pD,3))];    
    bar([mean(cD1,2),mean(pD1,2)]);
    
    hold on
    
    % Finding the number of groups and the number of bars in each group
    ngroups = size(mRes,1);
    nbars = size(mRes,2);

    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, mRes(:,i), sRes(:,i), 'k', 'linestyle', 'none');
    end


    legend({'Control','Papa'})
    set(gca,'xtick',[1:length(roiLabel)],'xticklabel',roiLabel)
    set(gca,'TickLabelInterpreter','none')
    ylabel('%BOLD signal change','FontSize',14)
    title(contrast_no,'FontSize',14)
    xticklabel_rotate([1:length(roiLabel)],90,roiLabel,'interpreter','none','fontsize',11)
    if save_figure
        print(fh,fullfile(resultDir,[contrast_names{iCon} '_bar_' datestr(datetime('now'),'yymmddTHHMMSS') '.pdf']),'-dpdf')
        savefig(fh,fullfile(resultDir,[contrast_names{iCon} '_bar_' datestr(datetime('now'),'yymmddTHHMMSS') '.fig']))
    end

end

%% Save data in tables; one per Contrast

fileID = fullfile(resultDir,sprintf('%s_PvsC_percSignalChange_allROI_%s.xlsx',task,datestr(datetime('now'),'yymmddTHHMMSS')));

for iCon = 1:length(contrast_names)
    
    pt = p(iCon,:)';
    cD0 = filloutliers(squeeze(cD(iCon,:,:)),'clip');
    pD0 = filloutliers(squeeze(pD(iCon,:,:)),'clip');
    cD1 = mean(cD0,2);
    pD1 = mean(pD0,2);  
    cD2 = std(cD0,[],2)/sqrt(size(cD,3));
    pD2 = std(pD0,[],2)/sqrt(size(pD,3)); 
    
    potato = [pD1 pD2 cD1 cD2 pt]; 
    
    T = [table(roiLabel','variablenames',{'ROI'}) array2table(potato,'variablenames',{'Papa_mean','Papa_SE','Con_mean','Con_SE','P'})];    
        
    % save contrast estimates
    
    writetable(T,fileID,'Sheet',contrast_names{iCon});    
    
end

%% Save data in tables for individual subjects; one per Contrast

fileID = fullfile(resultDir,sprintf('%s_PvsC_percSignalChange_allROI_indvSubj_%s.xlsx',task,datestr(datetime('now'),'yymmddTHHMMSS')));

for iCon = 1:length(contrast_names)
    
    contrast_no = contrast_names{iCon};
    voxelDataSaveFile = dir(fullfile(voxelDataSaveDir,['ROI_data_joinROIs_' signal_mode '_' contrast_no '_*']));
    load(fullfile(voxelDataSaveFile(end).folder,voxelDataSaveFile(end).name))
    
    cIdx = find([S.Group]==-0.5);
    pIdx = find([S.Group]==0.5);
    
    % Save control subjects' data to an excel file
    vD = zeros(length(subjID_list(cIdx)),length(roiLabel));
    for iROI = 1:length(roiLabel)                
        roiVoxelData = reshape(mean([voxelData_all(cIdx,iROI).expVoxelData]),[],1)-reshape(mean([voxelData_all(cIdx,iROI).conVoxelData]),[],1);
        vD(:,iROI) = roiVoxelData;        
    end
    
    T = [table(subjNo_list(cIdx),'variablenames',{'Subj'}) array2table(vD,'variablenames',roiLabel)];    
    writetable(T,fileID,'FileType','spreadsheet','Sheet',contrast_names{iCon},'Range','A1');  
    
    % Save papa subjects' data to an excel file
    vD = zeros(length(subjID_list(pIdx)),length(roiLabel));
    for iROI = 1:length(roiLabel)                
        roiVoxelData = reshape(mean([voxelData_all(pIdx,iROI).expVoxelData]),[],1)-reshape(mean([voxelData_all(pIdx,iROI).conVoxelData]),[],1);
        vD(:,iROI) = roiVoxelData;        
    end
    
    T = [table(subjNo_list(pIdx),'variablenames',{'Subj'}) array2table(vD,'variablenames',roiLabel)];    
    writetable(T,fileID,'FileType','spreadsheet','Sheet',contrast_names{iCon},'Range',[cast('A'+length(roiLabel)+2,'char') '1']);
        
end

disp('Job''s done, boss.')