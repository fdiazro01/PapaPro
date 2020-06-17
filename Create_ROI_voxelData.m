clear, clc


% Implement preprocessing from BrainDecoderToolbox v2
%% Parameters

task = 'video';
group = {'P','C'};
smoothing = 's4';
TR = 3; 
nRun = 2; % Number of runs per protocold
nSlicesPerEpoch = 10;

% Import subjects' data

% Import the data
[S subjNo_list subjID_list] = read_subject_data("data.xlsx",group);

varNames = fieldnames(S);
nVar = length(varNames);

% Eliminate subjects that are to be eliminated (for example, due to large movement artifacts)
if strcmp(task,'audio')
    rIdx = find(~[S.Audio]);
    subjNo_list(rIdx) = [];
    subjID_list(rIdx) = [];
    S(rIdx) = [];   
elseif strcmp(task,'video')
    rIdx = find(~[S.Video]);
    subjNo_list(rIdx) = [];
    subjID_list(rIdx) = [];
    S(rIdx) = [];   
end


% Index to where the control and papa group data are
grIdx = unique([S.Group]);
cIdx = find([S.Group]==grIdx(1));
pIdx = find([S.Group]==grIdx(2));
    
% Directories and contrasts

workDir = []; % Directory where SPM 1st level analysis data is saved
voxelDataSaveDir = []; % Directory to save the voxel data

stimDir = []; % Directory where the stimuli presentation order information is saved
stimFileName = 'stimOrder_run';

roiDir = []; % Directory with the saved ROI mask files
roiFilter = 'mat_*.mat';

joinROIs = {'AMY','STS','TP'}; % Join the individual lateral parts of these ROIs together

save_data = 1;

skip_ifprocessed = 0;

signal_mode = 'PercentSignalChange';
% signal_mode = 'Zscore';

%% Preload ROI voxel data

clear roiData

roiFiles = dir(fullfile(roiDir,roiFilter));
for iROI = 1:length(roiFiles)
    roiName = roiFiles(iROI).name;    
    roiData_all(iROI) = load(fullfile(roiDir,roiName));
    idx = strfind(roiName,'_');
    roiLabel{iROI} = roiName(idx(1)+1:idx(3)-1);
end

% Join L and R parts of the ROI
if ~isempty(joinROIs)
    for i=1:length(joinROIs)
        j = joinROIs{i};
        idx = find(contains(roiLabel,j));
        roiData_all(idx(1)).cXYZ = [roiData_all(idx(1)).cXYZ,roiData_all(idx(2)).cXYZ];    
        idx2del(i) = idx(2);
        roiLabel{idx(1)} = j;
    end
    roiData_all(idx2del) = [];
    roiLabel(idx2del) = [];
end

%% Load data and calculate correlation between two ROIs

contrast_list = {'SvsC','S1vsC1','S2vsC2'};

for iCon=1:length(contrast_list)
    
    % Contrast to use. Either SvsC, S1vsC1, or S2vsC2
    contrast2use = contrast_list{iCon};

    if strcmp(contrast2use,'SvsC')
        % Index of exp and control stimulus
        stimExpIdx = [1 3];
        stimConIdx = [2 4];
    elseif strcmp(contrast2use,'S1vsC1')
        stimExpIdx = [1];
        stimConIdx = [2];
    elseif strcmp(contrast2use,'S2vsC2')
        stimExpIdx = [3];
        stimConIdx = [4];
    else
        error('Wrong contrast.\n')
    end
    
    Xexp_all = zeros(length(subjID_list),length(roiLabel));
    Xcon_all = zeros(length(subjID_list),length(roiLabel));
    Xbsl_all = zeros(length(subjID_list),length(roiLabel));
    X_all = zeros(length(subjID_list),length(roiLabel));
    
    for iSubj = 1:length(subjID_list)

        % Loads SPM.mat
        subjName = [subjNo_list{iSubj}(1:4) '_' subjID_list{iSubj}];
        load(fullfile(workDir,subjName,'SPM.mat'))

        expOnsets = [];
        conOnsets = [];    
        for i=1:nRun
            stim = load(fullfile(stimDir,subjName,task,[stimFileName num2str(i) '.mat']));
            onsets = cell2mat([stim(:).onsets]);        
            expOnsets = [expOnsets onsets(:,stimExpIdx)]; %#ok<*AGROW>
            conOnsets = [conOnsets onsets(:,stimConIdx)];

        end    
        d = stim.durations{1}/TR; % Durations MUST be the same. If not, change here
        expOnsets = expOnsets/TR; % Change units from seconds to slices
        conOnsets = conOnsets/TR; % Change units from seconds to slices

        nStimPerRun = size(expOnsets,2)/nRun;

        % Load movement parameter data
        Mfile = dir(fullfile(stimDir,subjName,task,'Nifti\*\*\rp_*.txt'));
        M = [];
        for i=1:length(Mfile)
            M = [M;load(fullfile(Mfile(i).folder,Mfile(i).name))];
        end

        % Load ROI data  
        M1 = M;
        for iROI = 1:length(roiLabel)

            roiName = roiLabel{iROI};    

            M = M1;
            fprintf('Extracting data for subj %s - ROI:%s...\n',subjName,roiName)

            roiData = roiData_all(iROI);
            voxelData = spm_get_data(SPM.xY.P, roiData.cXYZ);            
            nSlicesPerRun = size(voxelData,1)/nRun;

            % Stimulus indices
            idx = zeros(size(voxelData,1),1);
            for i=1:nRun
                eO = expOnsets(:,(i-1)*nStimPerRun+1:nStimPerRun*i); %#ok<*AGROW>
                eO = (i-1)*nSlicesPerRun+eO(:);
                eO = repmat(eO,1,d)+repmat([1:d],length(eO),1);
                idx(eO(:)) = 1;

                cO = conOnsets(:,(i-1)*nStimPerRun+1:nStimPerRun*i); %#ok<*AGROW>
                cO = (i-1)*nSlicesPerRun+cO(:);
                cO = repmat(cO,1,d)+repmat([1:d],length(cO),1);
                idx(cO(:)) = -1;
            end            

            % Shift the data 2 volumes (6s) to account for hemodynamic
            % delay
            [voxelData, ind] = shift_sample(voxelData, 'Group', [ones(size(voxelData,1)/2,1);ones(size(voxelData,1)/2,1)+1], 'ShiftSize', 2);
            idx = idx(ind);
            M = M(ind,:);

            % Regressout movement parameters, DC removal and linear
            % detrend, reduce outliers and normalize voxel data                
            voxelData = regressout(voxelData,'Group',[ones(size(voxelData,1)/2,1);ones(size(voxelData,1)/2,1)+1],...
                                    'Regressor',M,'LinearDetrend','on');               
            voxelData_norm = normalize_sample(voxelData, 'Group', [ones(size(voxelData,1)/2,1);ones(size(voxelData,1)/2,1)+1],'Mode',signal_mode,'Baseline',idx==0);
            voxelData_norm = reduce_outlier(voxelData_norm,'Dimension',2); 

            expVoxelData = voxelData_norm(idx==1,:); % Voxel data for the experimental condition
            conVoxelData = voxelData_norm(idx==-1,:);  % Voxel data for the control condition
            bslVoxelData = voxelData_norm(idx==0,:);  % Voxel data for the baseline

            for i=1:nRun*length(stimExpIdx)                    
                expVoxelData1(i,:,:) = expVoxelData(1+(i-1)*nSlicesPerEpoch:i*nSlicesPerEpoch,:);
                conVoxelData1(i,:,:) = conVoxelData(1+(i-1)*nSlicesPerEpoch:i*nSlicesPerEpoch,:);
            end

            for i=1:nRun*length(stimExpIdx)-1                    
                bslVoxelData1(i,:,:) = bslVoxelData(1+(i-1)*nSlicesPerEpoch:i*nSlicesPerEpoch,:);
            end

            expVoxelData_all = mean(expVoxelData,1);
            conVoxelData_all = mean(conVoxelData,1);
            bslVoxelData_all = mean(bslVoxelData,1);

            expVoxelData_allVol = expVoxelData;
            conVoxelData_allVol = conVoxelData;
            bslVoxelData_allVol = bslVoxelData;

            expVoxelData = squeeze(mean(mean(expVoxelData1,1),3))';
            conVoxelData = squeeze(mean(mean(conVoxelData1,1),3))';
            bslVoxelData = squeeze(mean(mean(bslVoxelData1,1),3))';

            clear expVoxelData1 conVoxelData1 bslVoxelData1
              
            voxelData_all(iSubj,iROI).expVoxelData = expVoxelData; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).conVoxelData = conVoxelData; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).bslVoxelData = bslVoxelData; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).expVoxelData_all = expVoxelData_all; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).conVoxelData_all = conVoxelData_all; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).bslVoxelData_all = bslVoxelData_all; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).expVoxelData_allVol = expVoxelData_allVol; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).conVoxelData_allVol = conVoxelData_allVol; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).bslVoxelData_allVol = bslVoxelData_allVol; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).Subj = subjName; %#ok<*SAGROW>
            voxelData_all(iSubj,iROI).ROI = roiName; %#ok<*SAGROW>

        end

    end
    
    if save_data
        if ~isempty(joinROIs)           
            save(fullfile(voxelDataSaveDir,['ROI_data_joinROIs_' signal_mode '_' contrast2use '_' datestr(datetime('now'),'yymmddTHHMMSS') '.mat']),'voxelData_all','roiLabel','S','subjID_list','subjNo_list','joinROIs')
        else
            save(fullfile(voxelDataSaveDir,['ROI_data_' signal_mode '_' contrast2use '_' datestr(datetime('now'),'yymmddTHHMMSS') '.mat']),'voxelData_all','roiLabel','S','subjID_list','subjNo_list','joinROIs')
        end
        
    end

end

