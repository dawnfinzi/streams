%% Uses Kendrick's NCSNR map to calculate and save the noise ceiling
% Saves NC values for each subject, hemisphere and number of trials
% averaged

clear all

fs_base = [nsd_datalocation '/freesurfer/subj%02d'];

for subjix = 1:8 %eight subjects
    clear rh_nc lh_nc
    
    fs_dir = sprintf(fs_base, subjix);
    data_dir = sprintf('%s/ppdata/subj%02d/nativesurface/betas_fithrf_GLMdenoise_RR/',nsd_datalocation('betas'),subjix);  
    
    rh_nc = load_mgh([data_dir 'rh.ncsnr.mgh']);
    lh_nc = load_mgh([data_dir 'lh.ncsnr.mgh']);
    
    for n_trials = 1:3
        clear rh_nc_t lh_nc_t
        
        k = 1/n_trials;
        rh_nc_t = 100 * ((rh_nc.^2)./(rh_nc.^2+k));
        lh_nc_t = 100 * ((lh_nc.^2)./(lh_nc.^2+k));
        
        nsd_savemgz(rh_nc_t, sprintf([data_dir 'rh.nc_%dtrials.mgh'],n_trials),fs_dir)
        nsd_savemgz(lh_nc_t, sprintf([data_dir 'lh.nc_%dtrials.mgh'],n_trials),fs_dir)
    end
end