clear, clc

%% Parameters

voxelSize = 3; % in mm

fROI_box_side = 15; % \ Widths of the box ROI in mm
fROI_box_side_vox = fROI_box_side/3; % \ Widths of the box ROI in voxels

%% Directories and contrasts

workDir = []; % Directory where SPM data is saved. One SPM file is needed to obtain the transformation matrix, to change from MNI space to matrix space.
resultDir = []; % Directory in which to save results

if ~exist(resultDir,'dir'), mkdir(resultDir), end
% if ~exist(marsbarDir,'dir'), mkdir(marsbarDir), end

%% Import ROI data

% Loads a table that contains a row for each ROI. Each row has the
% following fields: hemisphere, x, y, and z (corresponding to the peak
% coordinate in MNI space based from previous research)
load("video_ROIpeaks_taken_from_Abraham2014.mat")

%% Anatomical ROI comparisons

contrast2use = 11;

bw = floor(fROI_box_side_vox/2);

dummyFile = dir(fullfile(workDir,'*',sprintf('con_00%2d.nii',contrast2use)));
dummyV = spm_vol(fullfile(dummyFile(1).folder,sprintf('con_00%2d.nii',contrast2use)));

for iROI = 1:size(ROIpeaksS1,1)
        
    dummy = cellstr(ROIpeaksS1.Hemisphere(iROI));
    roiName = [ROIpeaksS1.ROI{iROI} '_' dummy{:}];     
    
    % Peak activation voxel in the given ROI
    cXYZmm = [ROIpeaksS1.x(iROI) ROIpeaksS1.y(iROI) ROIpeaksS1.z(iROI)];
    pXYZmm = cXYZmm;    
    
    % Convert from MNI space to voxel matrix
    cXYZ = [cXYZmm(:,1) cXYZmm(:,2) cXYZmm(:,3) ones(size(cXYZmm,1),1)]*(inv(dummyV.mat))';
    cXYZ(:,4) = [];
    cXYZ = round(cXYZ); 
    pXYZ = cXYZ;
    
    % Find the voxels that are in a 5x5x5 box around the maximum
    % activation point REGARDLESS of whether they are in the given
    % ROI or not
    [x y z] = meshgrid([cXYZ(1)-bw:cXYZ(1)+bw],[cXYZ(2)-bw:cXYZ(2)+bw],[cXYZ(3)-bw:cXYZ(3)+bw]);
    cXYZ = [x(:)';y(:)';z(:)'];  
    
    roiLabel = sprintf('%s_%i_%i_%i',roiName,pXYZmm(1),pXYZmm(2),pXYZmm(3));
    
    save(fullfile(resultDir,['mat_' roiLabel '_box_' num2str(fROI_box_side_vox) '_vox.mat']),'cXYZ','pXYZ','pXYZmm')  
end



