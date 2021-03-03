#import packages
import argparse
import pandas as pd
import h5py
import numpy as np
import scipy as sp
import scipy.stats as stats
import scipy.io
import nibabel.freesurfer.mghformat as mgh
import scipy.io
import itertools 
import pickle
import sys

utils_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/code/streams/utils/'
sys.path.append(utils_dir)

from rsm_utils import get_flat_lower_tri, get_reliability_data 

data_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/'
local_data_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/'

def fast_pearson(x,y):
    #faster, vectorized version
    xz = x - x.mean(axis=0)
    yz = y - y.mean(axis=0)
    xzss = (xz * xz).sum(axis=0)
    yzss = (yz * yz).sum(axis=0)
    r = np.matmul(xz.transpose(), yz) / (np.sqrt(np.outer(xzss, yzss))+np.finfo(float).eps) #add machine prec. to avoid any divide by 0
    return np.maximum(np.minimum(r, 1.0), -1.0) #for precision issues

def main(subjid, hemi):
    n_repeats = 3

    #get trial ids and mask        
    all_ids = []
    max_session = np.zeros(len([subjid]))
    for sidx, sid in enumerate([subjid]):

        data = pd.read_csv(data_dir+'nsddata/ppdata/subj'+ sid +'/behav/responses.tsv', sep='\t')

        max_session[sidx] = np.max(np.array(data['SESSION'])) 

        all_ids.append(np.array(data['73KID']))

    which_reps = []
    for sidx, sid in enumerate([subjid]):
        vals, idx_start, count = np.unique(all_ids[sidx], return_counts=True,
                                        return_index=True)
        which_reps.append(vals[count == n_repeats])

    least_trials = min(which_reps, key=len)

    id_nums_3reps = []
    mask_3reps = []
    for sidx, sid in enumerate([subjid]):

        data = pd.read_csv(data_dir+'nsddata/ppdata/subj'+ sid +'/behav/responses.tsv', sep='\t')

        mask_3reps.append(np.isin(all_ids[sidx],which_reps[sidx]))
        id_nums_3reps.append(np.array(data['73KID'])[mask_3reps[sidx]])

    arr1inds = id_nums_3reps[sidx].argsort()
    
    rh_streams = []
    for sidx, sid in enumerate([subjid]):
        mgh_file = mgh.load(data_dir+'nsddata/freesurfer/subj'+ sid +'/label/rh.streams_shrink5.mgz')
        rh_streams.append(mgh_file.get_fdata()[:,0,0])
        
    stream_idx = np.where(rh_streams[0] != 0)
        
    for sidx, sid in enumerate([subjid]):
        
        mask = mask_3reps[sidx]
        sorted_betas = []

        #get all betas across all sessions
        for sess in range(1,int(max_session[sidx])+1):
            print(sess)

            if(sess < 10):
                idx = '0' + str(sess)
            else:
                idx = str(sess)

            raw_betas = h5py.File(local_data_dir+'freesurfer/subj'+sid+'/betas/'+ hemi +'.zscore_betas_session'+idx+'.hdf5','r')

            sess_betas = raw_betas['zscore_betas'][:][mask[(sess-1)*750:sess*750]]
            del raw_betas

            if(sess==1):
                betas_trimmed = sess_betas[:,stream_idx[0]]
            else:
                betas_trimmed = np.append(betas_trimmed,sess_betas[:,stream_idx[0]],axis=0)

            del sess_betas
        
    betas_by_repeat = []
    sorted_betas = betas_trimmed[arr1inds[::-1]]

    for r in range(n_repeats):
        betas_by_repeat.insert(r,sorted_betas[r::3])
    
    corr_struct_r1r2 = fast_pearson(betas_by_repeat[0],betas_by_repeat[1])
    
    #create dict for matlab
    cs = {}
    cs['matrix'] = corr_struct_r1r2[1:100,:]
    cs['idx'] = stream_idx

    #save out
    save_dir = '../../../local_data/processed'
    scipy.io.savemat(save_dir + '/r1r2corrs_mini.mat', cs)
    
if __name__ == "__main__":
    # Parse command line args
    parser = argparse.ArgumentParser()
    parser.add_argument("--subjid", type=str)
    parser.add_argument("--hemi", type=str)

    ARGS, _ = parser.parse_known_args()

    main(
        ARGS.subjid,
        ARGS.hemi,
    )
