"""
General utils for project
(credit for many of these to internal lab utils: https://github.com/neuroailab/)
"""

import numpy as np


def norm_image(x):
    return (x - np.min(x))/np.ptp(x)

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    exp = np.exp(x)
    return exp/ exp.sum(0) #sums over axis representing columns

def rsquared(predicted, actual):
    """The "rsquared" metric
    """
    a_mean = actual.mean()
    num = np.linalg.norm(actual - predicted)**2
    denom = np.linalg.norm(actual - a_mean)**2
    return 1 - num / denom

def featurewise_norm(data, fmean=None, fvar=None):
    """perform a whitening-like normalization operation on the data, feature-wise
       Assumes data = (K, M) matrix where K = number of stimuli and M = number of features
    """
    if fmean is None:
        fmean = data.mean(0)
    if fvar is None:
        fvar = data.std(0)
    data = data - fmean  #subtract the feature-wise mean of the data
    data = data / np.maximum(fvar, 1e-5)  #divide by the feature-wise std of the data
    return data, fmean, fvar