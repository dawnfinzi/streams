
subjid = 'subj05';
searchlight_fname = [nsd_datalocation('local') '/searchlight/' sprintf('searchix_radius16_fsaverage3_%s_hemi2.mat', subjid)];
load(searchlight_fname)

roivals = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/rh.streams.mgz',subjid));  % load in an existing file?

len_hemi = max(max(searchix));
roi_form = zeros([len_hemi 1]);

for ix = 1:len_hemi
    if roivals(ix) ~= 0 %needs to be in stream rois
        [a,b] = find(searchix == ix,1);
        if numel(a) ~= 0
            roi_form(ix,1) = a;
        end
    end
end
c = unique(roi_form);
ph = roi_form;
for t = 1:length(c)
    idx = find(ph == c(t));
    ph(idx) = t-1;
end
    


roivals = ph;    
roilabels=[];
cmap   = jet(256);
rng    = [0 max(ph)];
threshs = {[] [] []};

mgznames = {'corticalsulc' 'Kastner2015'};
crngs = {[0 28] [0 25]};
cmaps = {jet(256) jet(256)};
cvndefinerois;