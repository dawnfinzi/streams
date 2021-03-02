% load custom map
nc = load(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh_split_half.mat', subjid));
lh_nc = nc.mean';
nc = load(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/rh_split_half.mat', subjid));
rh_nc = nc.mean';

% background maps 
mgznames = {'corticalsulc' 'Kastner2015', {lh_nc rh_nc}};  % quantities of interest (1 x Q)
crngs = {[0 28] [0 25], [0.05 .5]};  % ranges for the quantities (1 x Q)
cmaps = {jet(256) jet(256), jet};  % colormaps for the quantities (1 x Q)
threshs = {0.5 0.5 [.05]};

% roi info
cmap   = jet(256);   % colormap for ROIs
rng    = [0 7];      % should be [0 N] where N is the max ROI index
roilabels = {'early' 'midventral' 'midlateral' 'midparietal' 'ventral' 'parietal' 'lateral' };  % 1 x N cell vector of strings

%% include previous rois drawn (or not)

% load in an existing file (both hemispheres)
roivals = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/?h.streams_shrink5.mgz',subjid));  % load in an existing file?

% load in rois for only one hemi (example: right)
%roivals = zeros(n_vertices,1);
%roivals((n_vertices/2)+1:end) = cvnloadmgz('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/fsaverage/label/rh.streams.mgz'); 

% blank slate
%roivals = [];

%% do it
cvndefinerois;
