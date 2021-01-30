clear all
close all

subjid = 'subj08';
searchlight_fname = [nsd_datalocation('local') '/searchlight/' sprintf('searchix_radius4_fsaverage3_%s_hemi2.mat', subjid)];
load(searchlight_fname)

roivals = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/rh.streams.mgz',subjid));  % load in an existing file?

len_hemi = size(roivals,1);
roi_form = zeros([len_hemi 1]);

count = 1;
for r = 1:size(searchix,1)
    row_r = searchix(r,:);
    roi = row_r(row_r~=0); %all searchlight idx for this roi
    
    streamvals = roivals(roi);
    percent_stream = sum(streamvals~=0)/size(streamvals,1);
    
    %if not overlapping with stream rois at all, get rid of it
    if percent_stream == 0 
        roi_form(roi) = 0;
    else
        roi_form(roi) = count;
        count = count+1;
    end
end

left = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/lh.radius4fsaverage3.mgz',subjid));  % load in an existing file?
roivals = [left; roi_form];
% % load in rois for only one hemi (example: right)
% n_vertices = length(roivals)*2;
% roivals = zeros(n_vertices,1);
% roivals((n_vertices/2)+1:end) = roi_form;
% 

%roivals = roi_form;    
roilabels=[];
cmap   = jet(256);
rng    = [0 max(roivals)];
threshs = {[] [] []};

mgznames = {'corticalsulc' 'Kastner2015'};
crngs = {[0 28] [0 25]};
cmaps = {jet(256) jet(256)};
cvndefinerois;