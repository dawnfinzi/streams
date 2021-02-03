%% Load and resave split-half reliability maps (created in python)

clear all

fs_base = [nsd_datalocation '/freesurfer/subj%02d'];
local_base = [nsd_datalocation('local') '/freesurfer/subj%02d'];

for subjix = 1:8 %eight subjects
    clear rh_sh lh_sh
    
    local_dir = sprintf(local_base, subjix);
    fs_dir = sprintf(fs_base, subjix);
    
    rh_sh = load([local_dir '/rh_split_half.mat']);
    rh_sh = rh_sh.mean';
    
    lh_sh = load([local_dir '/lh_split_half.mat']);
    lh_sh = lh_sh.mean';
    
    nsd_savemgz(rh_sh, [local_dir '/rh.avg_split_half.mgh'],fs_dir)
    nsd_savemgz(lh_sh, [local_dir '/lh.avg_split_half.mgh'],fs_dir)
    
end