"""
Z-score all betas per voxel and per session
"""

#import packages
import pandas as pd
import h5py
import numpy as np
import scipy as sp
import scipy.stats as stats
import scipy.io

#setup paths etc
data_dir = '../../../data/'
save_dir = '../../../local_data/freesurfer/'
subjid = ['01', '02', '03', '04','05', '06', '07', '08']
hemis = ['rh', 'lh']

## Calculate z-scored betas and save
for h, hemi in enumerate(hemis):

    for sidx, sid in enumerate(subjid):
    
        print(sid)
        
        data = pd.read_csv('../../../data/nsddata/ppdata/subj'+ sid +'/behav/responses.tsv', sep='\t')
        max_session = np.max(np.array(data['SESSION'])) 
        
        #get all betas across all sessions
        for sess in range(1,int(max_session)+1):
            
            if(sess < 10):
                idx = '0' + str(sess)
            else:
                idx = str(sess)

            raw_betas = h5py.File(data_dir+'nsddata_betas/ppdata/subj'+ sid +'/nativesurface/betas_fithrf_GLMdenoise_RR/'+hemi+'.betas_session'+idx+'.hdf5', 'r')
            betas = raw_betas['betas'][:]/300 
            
            z_betas = stats.zscore(betas) #zscore betas per voxel and session
            np.nan_to_num(z_betas,copy = False) #replace nans with 0, some voxels have 0 %signal change and we want to preserve that instead of turning to nans
            
            #save to local data folder
            h5f = h5py.File(save_dir+'subj'+sid+'/betas/'+hemi+'.zscore_betas_session'+idx+'.hdf5', 'w')
            h5f.create_dataset('zscore_betas', data=z_betas)
            h5f.close()
            
            #cleanup
            del raw_betas, betas, z_betas