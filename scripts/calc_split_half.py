  
"""
Calculate and save maps of split-half reliability for each voxel

Split-half reliability here is defined as the average correlation across the three combinations of pairs of repeats (i.e. presentation 1 correlated with presentation 2, presentation 2 with presentation 3, presentation 1 with presentation 3). The numbers of images presented 3 times varies slightly across participants.
"""

#import packages
import pandas as pd
import h5py
import numpy as np
import scipy as sp
import scipy.stats as stats
import scipy.io

data_dir = '../../../data/'
subjid = ['01', '02', '03', '04', '05', '06', '07', '08']
ROI_names = ['Unknown', 'Early', 'Midventral', 'Midlateral', 'Midparietal', 'Ventral', 'Lateral', 'Parietal']
n_repeats = 3
hemis = ['rh', 'lh']
#trial combos
t1 = [0, 0, 1]
t2 = [1, 2, 2]

## Find trials for which there were at least 3 repeats for each subjects
all_ids = []
max_session = np.zeros(len(subjid))
for sidx, sid in enumerate(subjid):
    data = pd.read_csv('../../../data/nsddata/ppdata/subj'+ sid +'/behav/responses.tsv', sep='\t')
    max_session[sidx] = np.max(np.array(data['SESSION'])) 
    all_ids.append(np.array(data['73KID']))

which_reps = []
for sidx, sid in enumerate(subjid):
    vals, idx_start, count = np.unique(all_ids[sidx], return_counts=True, return_index=True)
    which_reps.append(vals[count == n_repeats])    
least_trials = min(which_reps, key=len)

#use this to create a mask for the betas
mask_3reps = []
for sidx, sid in enumerate(subjid):
    data = pd.read_csv('../../../data/nsddata/ppdata/subj'+ sid +'/behav/responses.tsv', sep='\t')
    mask_3reps.append(np.isin(all_ids[sidx],which_reps[sidx]))


## Calculate and save out split-half for each hemisphere
for h, hemi in enumerate(hemis):
    for sidx, sid in enumerate(subjid):
        
        print(sid)
        mask = mask_3reps[sidx]
        
        #get all betas across all sessions
        for sess in range(1,int(max_session[sidx])+1):
            
            if(sess < 10):
                idx = '0' + str(sess)
            else:
                idx = str(sess)

            raw_betas = h5py.File(data_dir+'nsddata_betas/ppdata/subj'+ sid +'/nativesurface/betas_fithrf_GLMdenoise_RR/'+hemi+'.betas_session'+idx+'.hdf5', 'r')
            #betas = raw_betas['betas'][:]/300 #takes too much memory and conversion not necc for split half calcs
            
            sess_betas = raw_betas['betas'][:][mask[(sess-1)*750:sess*750]]
            del raw_betas

            if(sess==1):
                subj_betas = sess_betas
            else:
                subj_betas = np.concatenate((subj_betas, sess_betas))
                
            del sess_betas
        
        #sort betas into 1st, 2nd and 3rd presentations
        masked_ids = all_ids[sidx][mask]
        arr1inds = masked_ids.argsort()
        sorted_betas = subj_betas[arr1inds[::-1]]
        del subj_betas

        betas_by_repeat = []
        for r in range(n_repeats):
            betas_by_repeat.insert(r,sorted_betas[r::3])

        del sorted_betas

        n_vox = betas_by_repeat[0].shape[1]

        #calculate split-half reliability
        corrvals = np.zeros((n_vox,3))
        for vox in range(n_vox):
            for r in range(3):
                corrval = stats.pearsonr(betas_by_repeat[t1[r]][:,vox],
                                        betas_by_repeat[t2[r]][:,vox])[0]
                corrvals[vox, r] = corrval


        avg_corrvals = np.mean(corrvals, axis=1)
        sem_corrvals = stats.sem(corrvals, axis=1)

        #create dict for matlab
        split_half = {}
        split_half['mean'] = avg_corrvals
        split_half['sem'] = sem_corrvals

        #save out
        save_dir = '../../../local_data/freesurfer/subj' + sid 
        scipy.io.savemat(save_dir + '/'+hemi+'_split_half.mat', split_half)
        
        #cleanup
        del corrvals, avg_corrvals, sem_corrvals, split_half, betas_by_repeat
