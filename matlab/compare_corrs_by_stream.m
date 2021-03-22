clear all
close all
 
hh = 2;

subjids = {'subj01', 'subj02', 'subj03', 'subj04', ...
           'subj05', 'subj06', 'subj07', 'subj08'};
sids = {'01', '02', '03', '04', '05', '06', '07', '08'};       

for subjix = 1:length(subjids)

    subjid=subjids{subjix};
    sid = sids{subjix};

    parcelvals = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/freesurfer/%s/rh.tessellate_1000_trim10.mgz',subjid));  % load in an existing file?
    streamvals = cvnloadmgz(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/nsddata/freesurfer/%s/label/rh.streams.mgz',subjid));  % load in an existing file?

    load(sprintf('/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/processed/parcel_megas/%s_1000_tt10_5000imgs_rh_mega_matrix.mat',sid));
    mega_matrix=squeeze(matrix);

    %% determine stream roi for each parcel roi
    parcel_by_stream = [];

    for r = 1:length(mega_matrix)

        vals = streamvals(parcelvals == r); 
        vals = vals(vals~=0);

        %check that parcel has voxels only within a single stream ROI
        idx=find(all(~diff(vals)));
        if isempty(idx)
            parcel_by_stream(r) = 0;
        else
            parcel_by_stream(r) = vals(1,idx);
        end
    end

    lat_idx = find(parcel_by_stream == 7);
    lat_poss_perms = nchoosek(lat_idx,2);

    par_idx = find(parcel_by_stream == 6);
    par_poss_perms = nchoosek(par_idx,2);

    vent_idx = find(parcel_by_stream == 5);
    vent_poss_perms = nchoosek(vent_idx,2);

    l2p = combvec(lat_idx,par_idx)';
    l2v = combvec(lat_idx,vent_idx)';
    v2p = combvec(vent_idx, par_idx)';

    %%
    lat_corrs = [];
    for i = 1:length(lat_poss_perms)
        lat_corrs = [lat_corrs mega_matrix(lat_poss_perms(i,1),lat_poss_perms(i,2))];
    end
    par_corrs = [];
    for i = 1:length(par_poss_perms)
        par_corrs = [par_corrs mega_matrix(par_poss_perms(i,1),par_poss_perms(i,2))];
    end
    vent_corrs = [];
    for i = 1:length(vent_poss_perms)
        vent_corrs = [vent_corrs mega_matrix(vent_poss_perms(i,1),vent_poss_perms(i,2))];
    end

    %%
    l2p_corrs = [];
    for i = 1:length(l2p)
        l2p_corrs = [l2p_corrs mega_matrix(l2p(i,1),l2p(i,2))];
    end

    v2p_corrs = [];
    for i = 1:length(v2p)
        v2p_corrs = [v2p_corrs mega_matrix(v2p(i,1),v2p(i,2))];
    end

    l2v_corrs = [];
    for i = 1:length(l2v)
        l2v_corrs = [l2v_corrs mega_matrix(l2v(i,1),l2v(i,2))];
    end
    
    within_single(subjix,:) = [mean(vent_corrs), mean(lat_corrs), mean(par_corrs)];
    across(subjix,:) = [mean(v2p_corrs), mean(l2v_corrs), mean(l2p_corrs)];
    
    haak_fig(subjix,:) = [mean([lat_corrs vent_corrs]), mean(l2v_corrs), ...
             mean([vent_corrs par_corrs]), mean(v2p_corrs), ...
             mean([lat_corrs par_corrs]), mean(l2p_corrs)];
    
end

%% within vs across streams
v = nanmean(haak_fig, 1);
values = reshape(v,[2,3])';
e = nanstd(haak_fig, 1) ./ sqrt(sum(~isnan(haak_fig), 1));
err = reshape(e,[2,3])';

niceFig([.1 .1 .6 .6])
caption_y='RSM correlation';
niceGroupedBars(values,err,{'Ventral/Lateral','Ventral/Parietal', 'Lateral/Parietal'},...
                {'Within','Between'})
ylim([0 0.4]);
ylabel(caption_y)

niceSave([nsd_datalocation('head') '/results/figures/'], ['haak_fig'])
% plot_vals = [mean([lat_corrs vent_corrs]), mean(l2v_corrs), ...
%              mean([lat_corrs par_corrs]), mean(l2p_corrs), ...
%              mean([vent_corrs par_corrs]), mean(v2p_corrs)];
%mybar(plot_vals)


% within_single = [mean(vent_corrs), mean(lat_corrs), mean(par_corrs)];
% across = [mean(v2p_corrs), mean(l2v_corrs), mean(l2p_corrs)];
%mybar([within_single across])
