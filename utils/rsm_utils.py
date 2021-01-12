"""
Utils for working with RSMs
"""

import numpy as np
import scipy.io
import nibabel.freesurfer.mghformat as mgh
import pickle


def get_flat_lower_tri(x, diagonal=False):
    """
    Returns the flattened lower triangle of a provided matrix
    Inputs
        x (np.ndarray): 2D matrix to get triangle from
        diagonal (bool): if True, keeps the diagonal as part of lower triangle
    """
    k = 0 if diagonal else -1
    lower_tri = np.tril(x, k).T.ravel()
    return lower_tri[lower_tri != 0]

def get_ROI_data(subjects, hemi):
    """
    Returns indices for all stream ROIs for requested subjects
    Inputs
        subjects: list of subject ids
        hemi: which hemisphere to pull data for ('rh' or 'lh')
    """
    streams = []
    for sidx, sid in enumerate(subjects):
        mgh_file = mgh.load(data_dir+'nsddata/freesurfer/subj'+ sid +'/label/'+ hemi +'.streams.mgz')
        streams.append(mgh_file.get_fdata()[:,0,0])
    return streams

def get_reliability_data(subjects, hemi):
    """
    Returns voxel-level split half reliability data for requested subjects and hemi
    Inputs
        subjects: list of subject ids
        hemi: which hemisphere to pull data for ('rh' or 'lh')
    """
    reliability = []
    for sidx, sid in enumerate(subjects):

        sh_dir = local_data_dir + 'freesurfer/subj' + sid + '/'+ hemi +'_split_half.mat'
        sh = scipy.io.loadmat(sh_dir)

        reliability.append(sh['mean'])
    return reliability

def make_flat_rsms(subjects, ROIs, hemi, thresh, local_data_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/Dawn/NSD/local_data/', zscore=True)):
    """
    Returns flattened lower triangle of RSM
    (RSM is RSM for each requested subject and ROI on the 515 shared images)
    Inputs
        subjects: list of subject ids
        ROIs: which ROIs to include
        hemi: which hemisphere to pull data for ('rh' or 'lh')
        thresh: what cutoff split-half reliability threshold to use for voxel inclusion
        local_data_dir: the data directory to pull from
        zscore: whether to use z-scored betas or non z-scored betas
    """
    #get ROI data
    streams = get_ROI_data(subjects, hemi)

    #get voxel level split-half reliability data
    reliability = get_reliability_data(subjects, hemi)

    #get organized betas
    if zscore:
        with open(local_data_dir + 'processed/'+ hemi +'_betas_by_repeat_by_ROI_zscore.data', 'rb') as filehandle:
            # read the data as binary data stream
            betas_by_repeat_by_ROI = pickle.load(filehandle)
    else:
        with open(local_data_dir + 'processed/'+ hemi +'_betas_by_repeat_by_ROI.data', 'rb') as filehandle:
            # read the data as binary data stream
            betas_by_repeat_by_ROI = pickle.load(filehandle)
        
    #Replace voxels with split-half reliability < thresh with NaNs and then trim those from data structure
    sh_by_ROI = [[[] for j in range(len(ROIs)-1)] for i in range(len(subjid))]
    #organize
    for sidx, sid in enumerate(subjid):  
        for roi_idx in range(len(ROIs)-1):       
            sh_by_ROI[sidx][roi_idx]=reliability[sidx][:,streams[sidx] == roi_idx+1]
    #convert to nans
    for sidx, sid in enumerate(subjid):  
        for roi_idx in range(len(ROIs)-1): 
            for vox in range(len(sh_by_ROI[sidx][roi_idx][0])):
                if sh_by_ROI[sidx][roi_idx][0][vox] < thresh:
                    betas_by_repeat_by_ROI[sidx][roi_idx][0][:,vox]=np.nan
                    betas_by_repeat_by_ROI[sidx][roi_idx][1][:,vox]=np.nan
                    betas_by_repeat_by_ROI[sidx][roi_idx][2][:,vox]=np.nan    
    #trim out nans
    for sidx, sid in enumerate(subjid):   
        for roi_idx in range(len(ROIs)-1): 
            for r in range(n_repeats):
                temp = betas_by_repeat_by_ROI[sidx][roi_idx][r]
                trimmed = temp[:,~np.all(np.isnan(temp), axis=0)]

                betas_by_repeat_by_ROI[sidx][roi_idx][r] = trimmed
   
    #Create RSMS for all the ROIs, repeats and subjects
    tril_flat_shape = int((betas_by_repeat_by_ROI[0][0][0].shape[0]**2/2) - (betas_by_repeat_by_ROI[0][0][0].shape[0]/2))
    flat_rsm = np.zeros((len(subjid),len(ROI_names)-1, tril_flat_shape, n_repeats))
    rsm = np.zeros((len(subjid),len(ROI_names)-1,n_repeats,betas_by_repeat_by_ROI[0][0][0].shape[0],betas_by_repeat_by_ROI[0][0][0].shape[0]))

    for sidx, sid in enumerate(subjid):
        for roi_idx in range(len(ROIs)-1):
            for r in range(n_repeats):
                rsm[sidx,roi_idx,r,:,:] = np.corrcoef(betas_by_repeat_by_ROI[sidx][roi_idx][r])
                flat_rsm[sidx, roi_idx, :,r] = get_flat_lower_tri(rsm[sidx,roi_idx,r,:,:],diagonal=False)
                #lower = np.tril(rsm[sidx,roi_idx,r,:,:], -1).T.ravel() #only need lower triangle without diagonal
                #flat_rsm[sidx, roi_idx, :,r] = lower[lower != 0] #flatten into vectors
    
    return flat_rsm