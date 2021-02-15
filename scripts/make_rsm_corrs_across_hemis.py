"""
Given some subject and set of ROIs in both hemispheres, create RSMs across all (~10000) images and correlate RSMs across ROIs

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

utils_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/code/streams/utils/'
sys.path.append(utils_dir)

from rsm_utils import get_flat_lower_tri, get_reliability_data 

data_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/data/'
local_data_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/'


def main(subjid, roi_name, thresh=0.0):

    n_repeats = 3
    hemis = ['lh', 'rh']
    both_flat_rsm = []
    both_num_rois = []

    for h, hemi in enumerate(hemis):
        #get ROI data
        parcels = []
        mgh_file = mgh.load(local_data_dir+'freesurfer/subj'+ subjid +'/' + hemi + '.' + roi_name + '.mgz')
        parcels.append(mgh_file.get_fdata()[:,0,0])

        num_rois = int(np.max(parcels))

        #get split-half reliablity for each voxel
        reliability = get_reliability_data([subjid], hemi)

        sh_by_ROI = [[[] for j in range(num_rois)] for i in range(len([subjid]))]
        total_vox = np.zeros((len([subjid]), num_rois))

        for sidx, sid in enumerate([subjid]):  
            for roi_idx in range(num_rois):       
                sh_by_ROI[sidx][roi_idx]=reliability[sidx][:,parcels[sidx] == roi_idx+1]
                total_vox[sidx,roi_idx] = len(sh_by_ROI[sidx][roi_idx][0])

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

        #Replace voxels with split-half reliability < thresh with NaNs and then trim those from data structure
        #convert to nans
        for sidx, sid in enumerate([subjid]):  
            for roi_idx in range(num_rois): 
                for vox in range(len(sh_by_ROI[sidx][roi_idx][0])):
                    if sh_by_ROI[sidx][roi_idx][0][vox] < thresh:
                        betas_by_repeat_by_ROI[sidx][roi_idx][0][:,vox]=np.nan
                        betas_by_repeat_by_ROI[sidx][roi_idx][1][:,vox]=np.nan
                        betas_by_repeat_by_ROI[sidx][roi_idx][2][:,vox]=np.nan    
        #trim out nans
        for sidx, sid in enumerate([subjid]):   
            for roi_idx in range(num_rois): 
                for r in range(n_repeats):
                    temp = betas_by_repeat_by_ROI[sidx][roi_idx][r]
                    trimmed = temp[:,~np.all(np.isnan(temp), axis=0)]

                    betas_by_repeat_by_ROI[sidx][roi_idx][r] = trimmed
        
        #Create RSMS for all the ROIs, repeats and subjects
        tril_flat_shape = int((betas_by_repeat_by_ROI[0][0][0].shape[0]**2/2) - (betas_by_repeat_by_ROI[0][0][0].shape[0]/2))
        hemi_flat_rsm = np.zeros((num_rois, tril_flat_shape, n_repeats))

        sidx = 0 #currently doing one subject at a time
        for roi_idx in range(num_rois):
            for r in range(n_repeats):
                rsm = np.corrcoef(betas_by_repeat_by_ROI[sidx][roi_idx][r])
                hemi_flat_rsm[roi_idx, :, r] = get_flat_lower_tri(rsm,diagonal=False)

        both_flat_rsm.append(hemi_flat_rsm)
        both_num_rois.append(num_rois)

        del hemi_flat_rsm
        del betas_by_ROI
        del betas_by_repeat_by_ROI

    r1_trial_order = [0, 0, 1, 1, 2, 2]
    r2_trial_order = [1, 2, 0, 2, 0, 1]

    #make the mega matrix!
    total_combs = both_num_rois[0] + both_num_rois[1]
    mega_matrix = np.zeros((total_combs,total_combs))

    for idx1 in range(total_combs): #rows - i.e. model candidate

        if idx1 < both_num_rois[0]:
            hidx1 = 0 #left hemi
            roi_idx1 = idx1
        else:
            hidx1 = 1 #right hemi
            roi_idx1  = int(idx1-both_num_rois[0]) #roi
        
        split_half = np.zeros((3))
        split_half = [stats.pearsonr(both_flat_rsm[hidx1][roi_idx1,:,0],both_flat_rsm[hidx1][roi_idx1,:,1])[0],
                    stats.pearsonr(both_flat_rsm[hidx1][roi_idx1,:,0],both_flat_rsm[hidx1][roi_idx1,:,2])[0],
                    stats.pearsonr(both_flat_rsm[hidx1][roi_idx1,:,1],both_flat_rsm[hidx1][roi_idx1,:,2])[0]]
        NC_model = np.mean(split_half) * 100
        
        for idx2 in range(total_combs): #columns - i.e. target data

            if idx2 < both_num_rois[0]:
                hidx2 = 0 #left hemi
                roi_idx2 = idx2
            else:
                hidx2 = 1 #right hemi
                roi_idx2  = int(idx2-both_num_rois[0]) #roi
            
            split_half = np.zeros((3))
            split_half = [stats.pearsonr(both_flat_rsm[hidx2][roi_idx2,:,0],both_flat_rsm[hidx2][roi_idx2,:,1])[0],
                        stats.pearsonr(both_flat_rsm[hidx2][roi_idx2,:,0],both_flat_rsm[hidx2][roi_idx2,:,2])[0],
                        stats.pearsonr(both_flat_rsm[hidx2][roi_idx2,:,1],both_flat_rsm[hidx2][roi_idx2,:,2])[0]]
            NC_target = np.mean(split_half) * 100
            
            rsm_corr = np.zeros((6))
            for r in range(6):
                rsm_corr[r] = stats.pearsonr(both_flat_rsm[hidx1][roi_idx1,:, r1_trial_order[r]],
                                            both_flat_rsm[hidx2][roi_idx2,:, r2_trial_order[r]])[0]
            
            mega_matrix[idx1,idx2] = np.mean(rsm_corr) * np.sqrt(100/NC_model) * np.sqrt(100/NC_target)

    #save to local data folder
    save_file = local_data_dir + 'processed/' + subjid + '_'  + roi_name + '_both.data'

    with open(save_file, 'wb') as filehandle:
        # store the data as binary data stream
        pickle.dump([mega_matrix], filehandle)


if __name__ == "__main__":
    # Parse command line args
    parser = argparse.ArgumentParser()
    parser.add_argument("--subjid", type=str)
    parser.add_argument("--roi_name", type=str)
    parser.add_argument("--thresh", type=float, default=0.0)
    ARGS, _ = parser.parse_known_args()

    main(
        ARGS.subjid,
        ARGS.roi_name,
        thresh = ARGS.thresh,
    )
