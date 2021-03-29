"""
Given some subject, hemisphere and set of ROIs, create RSMs across all (~10000) images and correlate RSMs across ROIs

Saves a "mega matrix" that is the correlation matrix for the ROIs
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
from timeit import default_timer as timer
utils_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/code/streams/utils/'
sys.path.append(utils_dir)

from rsm_utils import get_flat_lower_tri, get_reliability_data 

data_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/'
local_data_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/'

def vcorrcoef(X,y):
    Xm = np.reshape(np.mean(X,axis=1),(X.shape[0],1))
    ym = np.mean(y)
    r_num = np.sum((X-Xm)*(y-ym),axis=1)
    r_den = np.sqrt(np.sum((X-Xm)**2,axis=1)*np.sum((y-ym)**2))
    r = r_num/r_den
    return r

def main(subjid, hemi, roi_name):
    
    print(subjid)
    start = timer() #start the clock

    n_repeats = 3

    #get ROI data
    parcels = []
    mgh_file = mgh.load(local_data_dir+'freesurfer/subj'+ subjid +'/' + hemi + '.' + roi_name + '.mgz')
    parcels.append(mgh_file.get_fdata()[:,0,0])

    num_rois = int(np.max(parcels))

    #get trial ids and mask
    expdesign = scipy.io.loadmat(data_dir+'nsddata/experiments/nsd/nsd_expdesign.mat')
    sharedix = expdesign['sharedix']        
    all_ids = []
    shared_mask = []
    shared_id_nums = []
    shared_3_reps = []
    shared_id_nums_3reps = []
    shared_mask_3reps = []
    max_session = np.zeros(len([subjid]))
    for sidx, sid in enumerate([subjid]):
        
        data = pd.read_csv(data_dir+'nsddata/ppdata/subj'+ sid +'/behav/responses.tsv', sep='\t')
        
        max_session[sidx] = np.max(np.array(data['SESSION'])) 
        
        all_ids.append(np.array(data['73KID']))

        shared_mask.append(np.isin(all_ids[sidx],sharedix))
        shared_id_nums.append(np.array(data['73KID'])[shared_mask[sidx]])

        vals, idx_start, count = np.unique(shared_id_nums[sidx], return_counts=True,
                                        return_index=True)
        shared_3_reps.append(vals[count == n_repeats])
        
    least_trials = min(shared_3_reps, key=len)

    for sidx, sid in enumerate([subjid]):
        
        data = pd.read_csv(data_dir+'nsddata/ppdata/subj'+ sid +'/behav/responses.tsv', sep='\t')
            
        shared_mask_3reps.append(np.isin(all_ids[sidx],least_trials))
        shared_id_nums_3reps.append(np.array(data['73KID'])[shared_mask_3reps[sidx]])

    arr1inds = shared_id_nums_3reps[sidx].argsort()

    #get and sort z-scored betas
    betas_by_ROI = [[] for j in range(num_rois)]

    for sidx, sid in enumerate([subjid]):
        
        mask = shared_mask_3reps[sidx]
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
    
    del betas_by_ROI #mems cleanup 
    del sorted_betas
    
    #Create RSMS for all the ROIs, repeats and subjects
    tril_flat_shape = int((betas_by_repeat_by_ROI[0][0][0].shape[0]**2/2) - (betas_by_repeat_by_ROI[0][0][0].shape[0]/2))
    flat_rsm0 = np.zeros((num_rois, tril_flat_shape))
    flat_rsm1 = np.zeros((num_rois, tril_flat_shape))
    flat_rsm2 = np.zeros((num_rois, tril_flat_shape))

    sidx = 0 #currently doing one subject at a time
    for roi_idx in range(num_rois):
        for r in range(n_repeats):
            rsm = np.corrcoef(betas_by_repeat_by_ROI[sidx][roi_idx][r])
            if r == 0:        
                flat_rsm0[roi_idx, :] = get_flat_lower_tri(rsm,diagonal=False)
            elif r == 1:
                flat_rsm1[roi_idx, :] = get_flat_lower_tri(rsm,diagonal=False)
            elif r == 2:
                flat_rsm2[roi_idx, :] = get_flat_lower_tri(rsm,diagonal=False)
            
    del betas_by_repeat_by_ROI #mems cleanup 
    
    # get NC for each ROI
    NC = np.zeros((num_rois))
    for ridx in range(num_rois):
        split_half = np.zeros((3))
        split_half = [stats.pearsonr(flat_rsm0[ridx,:],flat_rsm1[ridx,:])[0],
                    stats.pearsonr(flat_rsm0[ridx,:],flat_rsm2[ridx,:])[0],
                    stats.pearsonr(flat_rsm1[ridx,:],flat_rsm2[ridx,:])[0]]
        NC[ridx] = np.abs(np.mean(split_half) * 100) #no neg noise ceiling
    
    print('starting mega matrix')
    #make the mega matrix!
    mega_matrix = np.zeros((num_rois,num_rois))

    for roi_idx1 in range(num_rois): #rows - i.e. model candidate

        row = np.zeros((6,num_rois))
        for r in range(6): #loop through combos

            if r == 0:
                y = flat_rsm0[roi_idx1,:] # 1 x k
                X = flat_rsm1 # N x k
            elif r == 1:
                y = flat_rsm0[roi_idx1,:] # 1 x k
                X = flat_rsm2 # N x k
            elif r == 2:
                y = flat_rsm1[roi_idx1,:] # 1 x k
                X = flat_rsm0 # N x k
            elif r == 3:
                y = flat_rsm1[roi_idx1,:] # 1 x k
                X = flat_rsm2 # N x k
            elif r == 4:
                y = flat_rsm2[roi_idx1,:] # 1 x k
                X = flat_rsm0 # N x k
            elif r == 5:
                y = flat_rsm2[roi_idx1,:] # 1 x k
                X = flat_rsm1 # N x k

            row[r,:] = vcorrcoef(X,y)

        mega_matrix[roi_idx1,:] = np.mean(row, axis = 0) * np.sqrt(100/NC[roi_idx1]) * np.sqrt(100/NC)


    #save to local data folder
    save_file = local_data_dir + 'processed/' + subjid + '_' + hemi + '_' + roi_name + '_sharedimgs_compbyrow.data'

    with open(save_file, 'wb') as filehandle:
        # store the data as binary data stream
        pickle.dump([mega_matrix], filehandle)
    
    end = timer()
    print(end - start)

if __name__ == "__main__":
    # Parse command line args
    parser = argparse.ArgumentParser()
    parser.add_argument("--subjid", type=str)
    parser.add_argument("--hemi", type=str)
    parser.add_argument("--roi_name", type=str)
    ARGS, _ = parser.parse_known_args()

    main(
        ARGS.subjid,
        ARGS.hemi,
        ARGS.roi_name,
    )
