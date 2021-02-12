function trim_out_rois(subjix, subjid, roi_name, thresh)

% function trim_out_rois(subix, subjid, roi_name, thresh)
%
% <subjix> is int
% <subjid> is is corresponding string for subject
% <thresh> is threshold for split-half reliability to use (float)
% <roi_name> is the inner name for the roi (string). i.e. rh.roi_name.mgz
%
% Example)
%subjix=5;
%subjid='subj05';
%thresh = .2;
%roi_name = 'tessellate_300';
%trim_out_rois(subjix, subjid, roi_name, thresh)
%
% code to take parcellation ROIs and trim out ROIs with no above threshold voxels
%

%% set up directories
fs_base = [nsd_datalocation '/freesurfer/subj%02d'];
fs_dir = sprintf(fs_base, subjix);

local_base = [nsd_datalocation('local') '/freesurfer/subj%02d'];
local_dir = sprintf(local_base, subjix);

%% read in data
%ROIs
left = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh.%s.mgz',subjid, roi_name));  % load in an existing file?
right = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/rh.%s.mgz',subjid, roi_name));  % load in an existing file?

%plit half reliability by voxel 
rh_sh = load([local_dir '/rh_split_half.mat']);
rh_sh = rh_sh.mean';

lh_sh = load([local_dir '/lh_split_half.mat']);
lh_sh = lh_sh.mean';

%% Left hemi deletion

% get num voxels above thresh for each ROI
for r = 1:max(left)
    sh_by_ROI = lh_sh(left==r);
    num_vox_l(r) = sum(sh_by_ROI > thresh);
end

zero_indices = find(num_vox_l == 0); %ROIs with no remaining voxels

% remove ROIs with no above threshold voxels
for l = 1:length(left)
    if sum(left(l)==zero_indices)>0 %roi with zero vox
        left(l) = 0;
    end
end


%% Reft hemi deletion

% get num voxels above thresh for each ROI
for r = 1:max(right)
    sh_by_ROI = rh_sh(right==r);
    num_vox_r(r) = sum(sh_by_ROI > thresh);
end

zero_indices = find(num_vox_r == 0);%ROIs with no remaining voxels

% remove ROIs with no above threshold voxels
for l = 1:length(right)
    if sum(right(l)==zero_indices)>0 %roi with zero vox
        right(l) = 0;
    end
end

%% Save out the new ROIs
nsd_savemgz(left, [local_dir '/lh.' roi_name '_trim' num2str(thresh*100) '.mgz'],fs_dir)
nsd_savemgz(right,  [local_dir '/rh.' roi_name '_trim' num2str(thresh*100) '.mgz'],fs_dir)
