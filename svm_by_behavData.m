% Same as Clustering_ROI_maxVoxels_video_svm.m except that it uses SEM
% instead of CI as errorbars. 

%v3: Classifies other behavioural data as well
%v4: Reads behavioural data from mat file saved during spm processing,
%    rather than from excel sheet, to prevent confusion upon updating the 
%    excel sheet.
%v5: Final version
%v5.05: Removes ROIs from the Mirror Neuron Network

clear, clc

%% Parameters

task = 'video';
group = {'P','C'};
% group = {'P'};
% group = {'C'};
smoothing = 'noSmooth'; % Do not use smoothing for multivoxel analysis

rois2join = {'AMY','STS','TP','vmPFC','LFPC'};

% SVM parameters
numVoxel = 0; % Use all voxels in the fROI

% - svmtrain
%     - `-t` : Kernel type (0 => linear)
%     - `-c` : Cost (parameter C of C-SVC)
%     - `-r` : Coefficient in kernel function
%     - `-d` : Degree in Kernel function
%     - `-b` : Tain SVC model for probaility estimates (1 => true)
%     - `-q` : Silent mode
% - svmpredict
%     - `-b` : Predict for probaility estimates (1 => true)
%
trainOpt = '-t 0 -c 1 -r 0 -d 3 -b 0 -q';
testOpt = '-b 0';

%% Directories and contrasts

resultDir = []; % Directory to save results
voxelDataSaveDir = []; % Directory in which the voxel data was saved

if ~exist(resultDir,'dir'), mkdir(resultDir), end

contrast_names = {'SvsC'};

nBoot = 100; % Number of repetitions of bootstrapping

rng('default') % Set like this to make the results of the bootstrapping reproduceable

%% Use SVM to classify, using voxeldata from S-C contrast

clear mR_base
    
contrast2use = contrast_names{1};
voxelDataSaveFile = 'ROI_data_joinROIs_PercentSignalChange_';
voxelDataFiles = dir(fullfile(voxelDataSaveDir,[voxelDataSaveFile contrast2use '_*']));
load(fullfile(voxelDataFiles(end).folder,voxelDataFiles(end).name))

dropROI_label = []; % Remove these ROIs from the analysis.
dropROI_idx = find(contains(roiLabel,dropROI_label));

roiLabel(dropROI_idx) = [];
voxelData_all(:,dropROI_idx) = [];

roiLabelJoined = roiLabel;
% Join L and R parts of the ROI
if ~isempty(rois2join)
    idx2del = [];
    for i=1:length(rois2join)
        j = rois2join{i};
        idx = find(contains(roiLabel,j));
        if isempty(idx) continue, end
        idx2del(i) = idx(2);
        roiLabelJoined{idx(1)} = j;
    end
    roiLabelJoined(idx2del) = [];
end

alpha = 0.05/length(roiLabelJoined); % Threshold for significance, adjusted for number of comparisons (number of ROIs)

grIdx = unique([S.Group]);
pIdx1 = find([S.Group]==grIdx(2));
cIdx1 = find([S.Group]==grIdx(1));

if any(strcmp(group,'P'))
    if any(strcmp(group,'C'))
        % When comparing between groups, only use the Group variable
        idx2use = [cIdx1,pIdx1];        
        var2use = [2]; 
    else
        % If using papa subjects, use all the non-redundant variables
        idx2use = pIdx1;        
        var2use = [3 4 14 15 16:18];
    end
elseif any(strcmp(group,'C'))
    % If using control subjects, remove the variables unique to fathers
    var2use = [3 14 16:18];
    idx2use = cIdx1;
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

binVar = {'Age','GA','Nmean','FetalAttachment','DAS','HouseIncome',...
           'Pmean','WeeklyWorktime','T','Oxy'}; % Binarize these variables by cutting them in half in the median (less than median = 0; higher than median = 1)

% Loads Parental Attitude scores
load('ParentAttitude_1st.mat')
subjNo_list1 = subjNo_list(idx2use);
potato = PA.Subj;
PA.logScore(~contains(potato,subjNo_list1),:) = [];
PA.meanScore(~contains(potato,subjNo_list1),:) = [];
PA.RawResponse(~contains(potato,subjNo_list1),:) = [];
PA.Subj(~contains(potato,subjNo_list1)) = [];

load('Hormone_1st.mat')
potato = H.SubjID;
rmIdx = find(~ismember(potato,subjNo_list(idx2use)));
H(rmIdx,:) = [];
potato(rmIdx) = [];
R(:,1) = H.Testosterone1st;
R(:,2) = H.Oxytocin1st;

if any(~strcmp(potato,subjNo_list(idx2use))) error('Hormonal data and subject data not aligned.'), end

% If there is any NaN value in the T levels, fill it with the theoretical lowest. 
lowestVal = 0.03;
R(isnan(R(:,1)),1) = lowestVal;

% Find unstable oxytocin data in order to remove it
rmIdx = find(H.Usability==0);
R(rmIdx,2) = NaN;

variables = {'Testosterone','Oxytocin'};

k = 1;
data(:,j+k) = PA.logScore(:,1);
varNames1{j+k} = 'Pmean';
k = k+1;
data(:,j+k) = PA.logScore(:,2);
varNames1{j+k} = 'Nmean';
k = k+1;

if any(strcmp(group,'P')) & ~any(strcmp(group,'C'))
    % Load corrected Fetal Attachment values
    load('FetalAttach_1st.mat')
    subjNo_list1 = subjNo_list(idx2use);
    potato = FA.Subj;
    FA.RawResponse(~contains(potato,subjNo_list1),:) = [];
    FA.Subj(~contains(potato,subjNo_list1)) = [];

    % Load corrected DAS (Dyadic Adjustment Scale) values
    load('DAS_1st.mat')
    subjNo_list1 = subjNo_list(idx2use);
    potato = DAS.Subj;
    DAS.RawResponse(~contains(potato,subjNo_list1),:) = [];
    DAS.Subj(~contains(potato,subjNo_list1)) = [];

    data(:,j+k) = mean(FA.RawResponse,2);
    varNames1{j+k} = 'FetalAttachment';
    k = k+1;
    data(:,j+k) = mean(DAS.RawResponse,2);
    varNames1{j+k} = 'DAS';
	k = k+1;
   
end

data(:,j+k) = R(:,1);
varNames1{j+k} = 'T';
k = k+1;
data(:,j+k) = R(:,2);
varNames1{j+k} = 'Oxy';
k = k+1;


%%
for iVar=1:size(data,2)

    targetIdx = data(:,iVar); % Use Group as the classifier

    data2(:,iVar) = targetIdx; % Keep track of the actual values before binarization
    
    if any(strcmp(varNames1{iVar},binVar))
       md = nanmedian(targetIdx);
       targetIdx(targetIdx<md) = 0; 
       targetIdx(targetIdx>=md) = 1;
    end  
    
    if strcmp(varNames1{iVar},'Group')
       uniqueVal = unique(targetIdx);
       targetIdx(targetIdx==uniqueVal(1)) = 0;
       targetIdx(targetIdx==uniqueVal(2)) = 1;
    end 
    
    data1(:,iVar) = targetIdx; % Keep track of values after binarization

    for iROI = 1:length(roiLabelJoined)

        roiName = roiLabelJoined{iROI};
        fprintf('Extracting data for %s - ROI:%s...\n',contrast2use,roiName)
        idx = find(contains(roiLabel,roiName));
        if length(idx)==1
            voxelData = (reshape([voxelData_all(:,idx).expVoxelData_all_masked],[],length(subjID_list))-reshape([voxelData_all(:,idx).conVoxelData_all_masked],[],length(subjID_list)))';
        else
            voxelData1 = (reshape([voxelData_all(:,idx(1)).expVoxelData_all_masked],[],length(subjID_list))-reshape([voxelData_all(:,idx(1)).conVoxelData_all_masked],[],length(subjID_list)))';
            voxelData2 = (reshape([voxelData_all(:,idx(2)).expVoxelData_all_masked],[],length(subjID_list))-reshape([voxelData_all(:,idx(2)).conVoxelData_all_masked],[],length(subjID_list)))';        
            voxelData = [voxelData1,voxelData2];
        end
        
        [~,rmIdx] = find(isnan(voxelData));
        voxelData(:,unique(rmIdx)) = [];
        
        rmIdx = find(isnan(targetIdx));
        targetIdx(rmIdx) = [];
        voxelData(rmIdx,:) = [];
        
        nFolds = length(targetIdx);

        result.cv.predScore = zeros(2,2); % Initializes confusion matrix for this fold
        
        for n = 1:nFolds

            fprintf('Fold %d\t', n);

            trainInd = 1:nFolds;
            trainInd(trainInd==n) = [];
            testInd = n;

            xTrain = voxelData(trainInd, :);
            xTest  = voxelData(testInd, :);
            tTrain = targetIdx(trainInd, :);
            tTest  = targetIdx(testInd, :);

            % Voxel selection based on ANOVA f values
            if numVoxel
                fVals = anova_fvals(xTrain, tTrain);
                xTrain = select_top(xTrain, fVals, numVoxel);
                xTest = select_top(xTest, fVals, numVoxel);
            end

            % Balance for number of subjects using weights
            if length(unique(tTrain))==2
                
                n1 = length(find(tTrain==0));
                n2 = length(find(tTrain==1));

                if n1>n2*1.1
                   xTrain_P = xTrain(find(tTrain==1),:);
                   xTrain_C = xTrain(find(tTrain==0),:);

                   bootResult = zeros(nBoot,1);
                   predLabel_list = zeros(nBoot,1);
                   for iBoot = 1:nBoot
                       xTrain_C1 = xTrain_C(randperm(n1,n2),:);

                       xTrain1 = [xTrain_P;xTrain_C1];
                       tTrain1 = [ones(n2,1);zeros(n2,1)];

                       % Training and test
                       trainResult = svmtrain(tTrain1, xTrain1, trainOpt);
                       [predLabel, predAccuracy, dummy] = svmpredict(tTest, xTest, trainResult, testOpt);

                       bootResult(iBoot) = predLabel == tTest;
                       predLabel_list(iBoot) = predLabel;

                   end

                   result.cv.fold(n).predAccuracy = 100*(mean(bootResult));
                   predLabel = round(mean(predLabel_list));
                   for i=1:length(predLabel)
                        result.cv.predScore(tTest(i)+1,predLabel(i)+1) = 1+result.cv.predScore(tTest(i)+1,predLabel(i)+1);
                   end
                   clear bootResult

                elseif n2>n1*1.1
                   xTrain_P = xTrain(find(tTrain==1),:);
                   xTrain_C = xTrain(find(tTrain==0),:);

                   bootResult = zeros(nBoot,1);
                   predLabel_list = zeros(nBoot,1);
                   for iBoot = 1:nBoot
                       xTrain_P1 = xTrain_P(randperm(n2,n1),:);

                       xTrain1 = [xTrain_P1;xTrain_C];
                       tTrain1 = [ones(n1,1);zeros(n1,1)];

                       % Training and test
                       trainResult = svmtrain(tTrain1, xTrain1, trainOpt);
                       [predLabel, predAccuracy, dummy] = svmpredict(tTest, xTest, trainResult, testOpt);

                       bootResult(iBoot) = predLabel == tTest;
                       predLabel_list(iBoot) = predLabel;

                   end

                   result.cv.fold(n).predAccuracy = 100*(mean(bootResult));
                   predLabel = round(mean(predLabel_list));
                   for i=1:length(predLabel)
                        result.cv.predScore(tTest(i)+1,predLabel(i)+1) = 1+result.cv.predScore(tTest(i)+1,predLabel(i)+1);
                   end
                   
                   clear bootResult                   

                else

                    trainResult = svmtrain(tTrain, xTrain, trainOpt);
                    [predLabel, predAccuracy, dummy] = svmpredict(tTest, xTest, trainResult, testOpt);   
                    result.cv.fold(n).predAccuracy = predAccuracy(1);
                    for i=1:length(predLabel)
                        result.cv.predScore(tTest(i)+1,predLabel(i)+1) = 1+result.cv.predScore(tTest(i)+1,predLabel(i)+1);
                    end
                end
            elseif length(unique(tTrain))==3
                n1 = length(find(tTrain==-1));
                n2 = length(find(tTrain==0));
                n3 = length(find(tTrain==1));
                
                [c,idx] = sort([n1 n2 n3]);
                c = c(1);
                xTrain_A = xTrain(find(tTrain==-1),:); %#ok<*FNDSB>
                xTrain_B = xTrain(find(tTrain==0),:);
                xTrain_C = xTrain(find(tTrain==1),:);
                
                bootResult = zeros(nBoot,1);
                for iBoot = 1:nBoot
                    xTrain_A1 = xTrain_A(randperm(n1,c),:);
                    xTrain_B1 = xTrain_B(randperm(n2,c),:);
                    xTrain_C1 = xTrain_C(randperm(n3,c),:);
                    xTrain1 = [xTrain_A1;xTrain_B1;xTrain_C1];
                    tTrain1 = [-1*ones(c,1);zeros(c,1);ones(c,1)];
                    % Training and test
                    trainResult = svmtrain(tTrain1, xTrain1, trainOpt);
                    [predLabel, predAccuracy, dummy] = svmpredict(tTest, xTest, trainResult, testOpt);
                    bootResult(iBoot) = predLabel == tTest;
                end

                result.cv.fold(n).predAccuracy = 100*round(mean(bootResult));
                result.cv.fold(n).predAccuracy = 100*(mean(bootResult));
                clear bootResult
                
            end

        end
        fprintf('Mean prediction accuracy = %.3f\n', mean([result.cv.fold(:).predAccuracy]));

        mR_base(iROI,iVar).predAccuracy = [result.cv.fold(:).predAccuracy]; 
        mR_base(iROI,iVar).mean = mean([result.cv.fold(:).predAccuracy]); 
        mR_base(iROI,iVar).sem = std([result.cv.fold(:).predAccuracy])/sqrt(nFolds); 
        mR_base(iROI,iVar).numClass = length(unique(targetIdx));
        mR_base(iROI,iVar).predScore = result.cv.predScore;
        mR_base(iROI,iVar).ROIsize = size(voxelData,2);
        mR_base(iROI,iVar).classSizes = [length(find(targetIdx==-1)) length(find(targetIdx==0)) length(find(targetIdx==1))];

        [mR_base(iROI,iVar).h mR_base(iROI,iVar).p mR_base(iROI,iVar).ci mR_base(iROI,iVar).stats] = ttest([result.cv.fold(:).predAccuracy],100/length(unique(targetIdx)),'Alpha',alpha);
        
        clear result
    end

end

%% Plot results
roiLabel = roiLabelJoined;
save_figure = 1;
focus_on_list = 1:length(varNames1);
% focus_on_list = [];
numClass = 2;
predAccuracy = reshape([mR_base(:).mean],size(mR_base));
predAccuracySem = reshape([mR_base(:).sem],size(mR_base));
predAccuracyPvalue = reshape([mR_base(:).p],size(mR_base));

fh = figure('PaperType', 'a4', 'PaperOrientation', 'landscape', ...
            'PaperUnit', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

ha = axes;
set(fh, 'CurrentAxes', ha);

hold on;
box off;

bar(ha, predAccuracy, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'BarWidth', 0.5, ...
    'ShowBaseLine', 'off');

ngroups = length(roiLabel);
nbars = size(predAccuracy,2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, predAccuracy(:,i), predAccuracySem(:,i), 'k', 'linestyle', 'none');
end

% Change level line
plot(ha, [0, length(roiLabel) + 1], [100 / numClass, 100 / numClass], ...
     'k--', ...
     'LineWidth', 1.0);

set(ha, 'LineWidth', 1.0, ...
        'Layer', 'top', ...
        'Clipping', 'off');
    
% Significance marks    
[sg_r sg_c] = ind2sub(size(mR_base),find([mR_base(:).h]));

for i=1:length(sg_r)
    text(sg_r(i)+(sg_c(i)-2)*0.25,predAccuracySem(sg_r(i),sg_c(i))+predAccuracy(sg_r(i),sg_c(i))+1,'*','FontSize',18,'Color','red')
end

title(['Classification using S-C contrast brain data (' group{:} ' group) - p<' num2str(alpha)],'interpreter','none');

legend(varNames1)
xlabel('ROI');
set(ha, 'XLim', [0 length(roiLabel)+1], 'XTick', [ 1:length(roiLabel) ], ...
        'XTickLabel', roiLabel,'ticklabelinterpreter','none','fontsize',14);
% set(xl,'interpreter','none')
    
ylabel('Prediction accuracy (%)');
set(ha, 'YLim', [0, 100]);

% Save name = task_base_group_contrast_errorbar_time
saveFileName = sprintf('%s_svm_results_(bootstrap)_%s_%s_sem_%s',task,contrast2use,group{:},datestr(datetime('now'),'yymmddTHHMMSS'));

if save_figure
   print(fh,fullfile(resultDir,[saveFileName '.pdf']),'-dpdf')
   savefig(fh,fullfile(resultDir,saveFileName))
   T = array2table(predAccuracy,'VariableNames',varNames1','RowNames',roiLabel);
   T1 = array2table(predAccuracySem,'VariableNames',varNames1','RowNames',roiLabel);
   T2 = array2table(predAccuracyPvalue,'VariableNames',varNames1','RowNames',roiLabel);
   writetable(T,fullfile(resultDir,saveFileName),'FileType','spreadsheet','WriteRowNames',true,'Sheet','MeanPredAcc');
   writetable(T1,fullfile(resultDir,saveFileName),'FileType','spreadsheet','WriteRowNames',true,'Sheet','SEMPredAcc');
   writetable(T2,fullfile(resultDir,saveFileName),'FileType','spreadsheet','WriteRowNames',true,'Sheet','PredAccPvalue');
end

if ~isempty(focus_on_list)
    
    for iFocus=1:length(focus_on_list)
        focus_on = focus_on_list(iFocus);
        numClass = mR_base(1,focus_on).numClass;
        predAccuracy = reshape([mR_base(:,focus_on).mean],size(mR_base,1),[]);
        predAccuracySem = reshape([mR_base(:,focus_on).sem],size(mR_base,1),[]);

        fh = figure('PaperType', 'a4', 'PaperOrientation', 'landscape', ...
                    'PaperUnit', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

        ha = axes;
        set(fh, 'CurrentAxes', ha);

        hold on;
        box off;

        bar(ha, predAccuracy, ...
            'LineStyle', 'none', ...
            'LineWidth', 1.0, ...
            'BarWidth', 0.5, ...
            'ShowBaseLine', 'off');
        errorbar(ha, predAccuracy, predAccuracySem, ...
                 'Color', 'k', ...
                 'LineStyle', 'none', ...
                 'LineWidth', 1.0, ...
                 'Marker', 'none');    

        % Change level line
        plot(ha, [0, length(roiLabel) + 1], [100 / numClass, 100 / numClass], ...
             'k--', ...
             'LineWidth', 1.0);

        set(ha, 'LineWidth', 1.0, ...
                'Layer', 'top', ...
                'Clipping', 'off');

        % Significance marks    
        sidx = find([mR_base(:,focus_on).p] < 0.05);
        for i=1:length(sidx)
            if mR_base(sidx(i),focus_on).h
                text(sidx(i)-0.25,predAccuracy(sidx(i))+predAccuracySem(sidx(i))+2,['* p=' num2str(mR_base(sidx(i),focus_on).p)],'FontSize',10,'Color','red')
            else
                text(sidx(i)-0.25,predAccuracy(sidx(i))+predAccuracySem(sidx(i))+2,['p=' num2str(mR_base(sidx(i),focus_on).p)],'FontSize',10,'Color','black')
            end
        end

        title([contrast_names{:} ' - svm classification - nVox: ' num2str(numVoxel) ' - ' group{:} ' group - ' varNames1{focus_on} '  - p<' num2str(alpha)...
               ' - class lengths: ' sprintf('[%s]',num2str(mR_base(1,iFocus).classSizes))],'interpreter','none');

        legend(varNames1{focus_on})
        xlabel('ROI');
        set(ha, 'XLim', [0 length(roiLabel)+1], 'XTick', [ 1:length(roiLabel) ], ...
                'XTickLabel', roiLabel,'ticklabelinterpreter','none','fontsize',14);
        
        ylabel('Prediction accuracy (%)');
        set(ha, 'YLim', [0, 100]);
        
        saveFileName = sprintf('%s_svm_results_(bootstrap)_%s_%s_%s_sem_%s',task,group{:},contrast2use,varNames1{focus_on},datestr(datetime('now'),'yymmddTHHMMSS'));

        if save_figure
           print(fh,fullfile(resultDir,[saveFileName '.pdf']),'-dpdf')
           savefig(fh,fullfile(resultDir,saveFileName))
        end
    end

end

disp('Job''s done.')