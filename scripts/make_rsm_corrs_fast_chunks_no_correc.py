"""
Given some subject, hemisphere and set of ROIs, create RSMs across all (~10000) images and correlate RSMs across ROIs

Saves a chunk of the "mega matrix" that is the correlation matrix for the ROIs
"""

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

def main(subjid, hemi, roi_name, min_idx, max_idx):
    
    chunks = list(range(min_idx,max_idx))

    print(chunks)
    print(subjid)

    n_repeats = 3
    convert_to_int = 10000

    #get ROI data
    parcels = []
    mgh_file = mgh.load(local_data_dir+'freesurfer/subj'+ subjid +'/' + hemi + '.' + roi_name + '.mgz')
    parcels.append(mgh_file.get_fdata()[:,0,0])

    num_rois = int(np.max(parcels))

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

    #get and sort z-scored betas
    betas_by_ROI = [[] for j in range(num_rois)]

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
                for roi_idx in range(num_rois):
                    betas_by_ROI[roi_idx] = sess_betas[:,parcels[sidx] == roi_idx+1]
            else:
                for roi_idx in range(num_rois):
                    betas_by_ROI[roi_idx] = np.append(betas_by_ROI[roi_idx],sess_betas[:,parcels[sidx] == roi_idx+1],axis=0)
            
            del sess_betas
    betas_by_repeat_by_ROI = [[[] for j in range(num_rois)] for i in range(len([subjid]))]
    for sidx, sid in enumerate([subjid]):
        for roi_idx in range(num_rois):  
            
            sorted_betas = betas_by_ROI[roi_idx][arr1inds[::-1]]
            
            for r in range(n_repeats):
                betas_by_repeat_by_ROI[sidx][roi_idx].insert(r,sorted_betas[r::3])
    
    #Create RSMS for all the ROIs, repeats and subjects (moved to within mega matrix loop for mem efficiency)
    tril_flat_shape = int((betas_by_repeat_by_ROI[0][0][0].shape[0]**2/2) - (betas_by_repeat_by_ROI[0][0][0].shape[0]/2))
    sidx = 0 #currently doing one subject at a time

    r1_trial_order = [0, 0, 1, 1, 2, 2]
    r2_trial_order = [1, 2, 0, 2, 0, 1]

    print('starting mega matrix')
    #make the mega matrix!
    mega_matrix = np.zeros((len(chunks),num_rois))
    m_idx = 0 #allow for chunks not starting at 0

    for cidx, roi_idx1 in enumerate(chunks): #rows - i.e. model candidate
        
        flat_rsm_roi1 = np.zeros((tril_flat_shape, n_repeats))
        for r in range(n_repeats):
            flip = betas_by_repeat_by_ROI[sidx][roi_idx1][r].T
            rsm = fast_pearson(flip,flip)
            flat_rsm_roi1[:, r] = get_flat_lower_tri(rsm,diagonal=False)
        
        flat_rsm_roi2 = np.zeros((tril_flat_shape, n_repeats)) #initialize once then rewrite for time
        rsm_corr = np.zeros((6))
        for roi_idx2 in range(num_rois): #columns - i.e. target data
            
            if betas_by_repeat_by_ROI[sidx][roi_idx2][0].size == 0: #no betas
                print(roi_idx2)
                mega_matrix[m_idx, roi_idx2] = 0
            else:
                for r in range(n_repeats):
                    flip = betas_by_repeat_by_ROI[sidx][roi_idx2][r].T
                    rsm = fast_pearson(flip,flip)
                    flat_rsm_roi2[:, r] = get_flat_lower_tri(rsm,diagonal=False)
                
                for r in range(6):
                    rsm_corr[r] = fast_pearson(flat_rsm_roi1[:, r1_trial_order[r]],
                                                flat_rsm_roi2[:, r2_trial_order[r]])[0][0]
                
                mega_matrix[m_idx,roi_idx2] = np.mean(rsm_corr) 
        m_idx += 1

    #save to local data folder
    save_file = local_data_dir + 'processed/' + subjid + '/' + hemi + '_' + roi_name + '_' + str(min_idx) + 'to' + str(max_idx) + '_nocorrection.data'

    with open(save_file, 'wb') as filehandle:
        # store the data as binary data stream
        pickle.dump([mega_matrix], filehandle)


if __name__ == "__main__":
    # Parse command line args
    parser = argparse.ArgumentParser()
    parser.add_argument("--subjid", type=str)
    parser.add_argument("--hemi", type=str)
    parser.add_argument("--roi_name", type=str)
    parser.add_argument("--min", type=int)
    parser.add_argument("--max", type=int)
    ARGS, _ = parser.parse_known_args()

    main(
        ARGS.subjid,
        ARGS.hemi,
        ARGS.roi_name,
        ARGS.min,
        ARGS.max,
    )
