%% 
% This script takes the extracted ROI voxel data and shows the average 
% across subjects grouped by behavioural variables, comparing whether 
% there are significant differences.

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


%% Functional ROI comparisons based on grouping from binarized variables

save_figure = 1;

clear h p ci stats
 
for iCon = 1:length(contrast_list)
    
    % Loads SPM.mat
    contrast_no = contrast_names{iCon};
    voxelDataSaveFile = dir(fullfile(voxelDataSaveDir,['ROI_data_joinROIs_' signal_mode '_' contrast_no '_*']));
    load(fullfile(voxelDataSaveFile(end).folder,voxelDataSaveFile(end).name))
    
    grIdx = unique([S.Group]);
    cIdx = find([S.Group]==grIdx(1));
    pIdx = find([S.Group]==grIdx(2));
    
    % Depending on the group, select only some variables
    if any(strcmp(group,'P'))
        if any(strcmp(group,'C'))
            idx2use = [cIdx,pIdx];
            var2use = [2 3 14 15 16:18];
        else
            idx2use = pIdx;
            var2use = [3 4 14 15 16:18];
        end
    elseif any(strcmp(group,'C'))
        var2use = [3 14 16:18];
        idx2use = cIdx;
    else
        error('Wrong group.')
    end

    S = S(idx2use);
    voxelData_all = voxelData_all(idx2use,:);
    subjID_list = subjID_list(idx2use);

    varNames = fieldnames(S);
    varNames1 = varNames(var2use);

    % Behavioral data
    data = zeros(length(S),length(var2use));
    for i=1:length(S)
        for j=1:length(var2use)
            data(i,j) = getfield(S,{i},varNames1{j});
        end
    end

    binVar = {'PE_Pets','PE_Child','RE_Child'}; % Binarize these variables by transforming one of the values into another
    binVar2 = {'Pmean','WeeklyWorktime','Age','GA','Nmean','Pbias','FetalAttachment','DAS','HouseIncome'}; % Binarize these variables by cutting them in half in the median (less than median = 0; higher than median = 1)
    
    % Loads updated Parental Attitude scores
    load('ParentAttitude_1st.mat')
    subjNo_list1 = subjNo_list(idx2use);
    temp = PA.Subj;
    PA.logScore(~contains(temp,subjNo_list1),:) = [];
    PA.meanScore(~contains(temp,subjNo_list1),:) = [];
    PA.RawResponse(~contains(temp,subjNo_list1),:) = [];
    PA.Subj(~contains(temp,subjNo_list1)) = [];

    data(:,j+1) = PA.logScore(:,1);
    varNames1{j+1} = 'Pmean';
    data(:,j+2) = PA.logScore(:,2);
    varNames1{j+2} = 'Nmean';
    data(:,j+3) = PA.logScore(:,1)./(PA.logScore(:,1)+PA.logScore(:,2));
    varNames1{j+3} = 'Pbias';

    if any(strcmp(group,'P')) & ~any(strcmp(group,'C'))
        % Load corrected Fetal Attachment values
        load('FetalAttach_1st.mat')
        subjNo_list1 = subjNo_list(idx2use);
        temp = FA.Subj;
        FA.RawResponse(~contains(temp,subjNo_list1),:) = [];
        FA.Subj(~contains(temp,subjNo_list1)) = [];

        % Load corrected DAS (Dyadic Adjustment Scale) values
        load('DAS_1st.mat')
        subjNo_list1 = subjNo_list(idx2use);
        temp = DAS.Subj;
        DAS.RawResponse(~contains(temp,subjNo_list1),:) = [];
        DAS.Subj(~contains(temp,subjNo_list1)) = [];

        data(:,j+4) = mean(FA.RawResponse,2);
        varNames1{j+4} = 'FetalAttachment';
        data(:,j+5) = mean(DAS.RawResponse,2);
        varNames1{j+5} = 'DAS';

    end
    
    nCol = 4;
    nRow = ceil(length(varNames1)/nCol);
    nFig = nCol*nRow;

    fh = figure('PaperType','a4','PaperOrientation','landscape',...
                        'PaperUnit','normalized','PaperPosition',[0 0 1 1]);
    ha = tight_subplot(nRow,nCol,[0.05 0.01],[0.1 0.1],[0.05 0.01]);                

    
    for iVar=1:length(varNames1)     
        
        clear cD pD
        
        targetIdx = data(:,iVar); % Use Group as the classifier

        if any(strcmp(varNames1{iVar},binVar))
           targetIdx(targetIdx==2) = 1; 
        end  
        if any(strcmp(varNames1{iVar},binVar2))
           md = median(targetIdx);
           targetIdx(targetIdx<md) = 0; 
           targetIdx(targetIdx>=md) = 1;
        end          

        cIdx = find(targetIdx==0);
        pIdx = find(targetIdx==1);

        for iROI = 1:length(roiLabel)

            fprintf('Extracting data for contrast %s - ROI:%s...\n',contrast_no,roiLabel{iROI})

            roiVoxelData = reshape(mean([voxelData_all(:,iROI).expVoxelData]),[],1)-reshape(mean([voxelData_all(:,iROI).conVoxelData]),[],1);
            cG = roiVoxelData(cIdx);
            pG = roiVoxelData(pIdx);

            [h(iCon,iVar,iROI),p(iCon,iVar,iROI),ci(iCon,iVar,iROI,:),stats(iCon,iVar,iROI)] = ttest2(pG,cG,'Alpha',0.05/length(roiLabel)); % Correct Alpha due to number of comparisons              

            cD(iCon,iROI,:) = cG;
            pD(iCon,iROI,:) = pG;

        end

        cD1 = squeeze(cD(iCon,:,:));
        pD1 = squeeze(pD(iCon,:,:));        
        
        cD1 = filloutliers(cD1,'clip');
        pD1 = filloutliers(pD1,'clip');
        
        mRes = [mean(cD1,2),mean(pD1,2)];
        sRes = [std(cD1,[],2)/sqrt(size(cD,3)),std(pD1,[],2)/sqrt(size(pD,3))];    
        
        bar(ha(iVar),[mean(cD1,2),mean(pD1,2)]);

        hold(ha(iVar),'on')

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
            errorbar(ha(iVar),x, mRes(:,i), sRes(:,i), 'k', 'linestyle', 'none');
        end

        legend({'<Median','>=Median'})
        
       
        title(ha(iVar),[varNames1{iVar} ' - ' contrast_no],'FontSize',14)
        set(ha(iVar),'ylim',[-0.2 0.5])
        if mod(iVar,nCol) == 1
            ylabel(ha(iVar),'%BOLD signal change','FontSize',14)
        else
            set(ha(iVar),'yticklabel',[])
%             ylabel('',14)
        end
                
        if iVar>nCol*(nRow-1)
            set(ha(iVar),'xtick',[1:length(roiLabel)],'xticklabel',roiLabel,'xticklabelrotation',90)
            set(ha(iVar),'TickLabelInterpreter','none','fontsize',12)
%             xticklabel_rotate([1:length(roiLabel)],90,roiLabel,'interpreter','none','fontsize',11)
        else
            set(ha(iVar),'xtick',[1:length(roiLabel)],'xticklabel',[],'fontsize',12)
        end
        
        p1 = squeeze(p(iCon,iVar,:));
        h1 = squeeze(h(iCon,iVar,:));
        sidx = find(p1 < 0.05);
        for i=1:length(sidx)
            if h1(sidx(i))
                text(ha(iVar),sidx(i)-0.25,max(mRes(sidx(i),:)+sRes(sidx(i),:)+0.02),['* p=' num2str(p1(sidx(i)))],'FontSize',10,'Color','red')
            else
                text(ha(iVar),sidx(i)-0.25,max(mRes(sidx(i),:)+sRes(sidx(i),:)+0.02),['p=' num2str(p1(sidx(i)))],'FontSize',10,'Color','black')
            end
        end
        
            
        
        saveFileName = sprintf('%s_%s_bar_%s',contrast_names{iCon},varNames1{iVar},datestr(datetime('now'),'yymmddTHHMMSS'));
        if save_figure
            print(fh,fullfile(resultDir,[saveFileName '.pdf']),'-dpdf')
            savefig(fh,fullfile(resultDir,[saveFileName '.fig']))
        end

    end

end
