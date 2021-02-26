function rois_tthresh(hemi, subjid, t_thresh)
% edited version of sonia's rois_tthresh to work particularly for stream
% rois
    % set FS directory
    fsdir = sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s', subjid); 

    % load custom map
    nc = load(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/%s_split_half.mat', subjid, hemi));
    thresh_by = nc.mean';

    % roi files
    input_file = sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/%s.streams.mgz',subjid, hemi);
    output_file = sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/%s.streams_shrink%s.mgz',subjid, hemi, num2str(100*t_thresh));

    % threshold rois
    roi0 = cvnloadmgz(input_file);
    for roi_num=1:max(roi0(:))% go through all of the ROIs in this file
        roi0(roi0==roi_num) = (thresh_by(roi0==roi_num)>t_thresh) * roi_num;
    end

    nsd_savemgz(roi0(:),output_file,fsdir);

    fprintf(['saving thresholded file in: ' output_file '...\n']);
end

