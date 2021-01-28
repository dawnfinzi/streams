close all
clear all

% set FS directory
setenv('SUBJECTS_DIR','/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer');

%% setup info for cvndefinerois
subjid = 'fsaverage';   % which subject

% load custom map
nc = load(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/results/%s/b3nc.mat', subjid));
n_vertices = length(nc.vals);
% divide into hemis for plotting
lh_nc = nc.vals(1:(n_vertices/2));
rh_nc = nc.vals((n_vertices/2)+1:end);

% background maps 
mgznames = {'corticalsulc' 'Kastner2015'}; {lh_nc rh_nc}};  % quantities of interest (1 x Q)
crngs = {[0 28] [0 25]}; [0 75]};  % ranges for the quantities (1 x Q)
cmaps = {jet(256) jet(256)}; jet};  % colormaps for the quantities (1 x Q)
threshs = {0.5 0.5 [10]};

% roi info
cmap   = jet(256);   % colormap for ROIs
rng    = [0 7];      % should be [0 N] where N is the max ROI index
roilabels = {'early' 'midventral' 'midlateral' 'midparietal' 'ventral' 'parietal' 'lateral' };  % 1 x N cell vector of strings

%% include previous rois drawn (or not)

% load in an existing file (both hemispheres)
roivals = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/?h.streams.mgz',subjid));  % load in an existing file?

% load in rois for only one hemi (example: right)
%roivals = zeros(n_vertices,1);
%roivals((n_vertices/2)+1:end) = cvnloadmgz('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/fsaverage/label/rh.streams.mgz'); 

% blank slate
%roivals = [];

%% do it
cvndefinerois;
