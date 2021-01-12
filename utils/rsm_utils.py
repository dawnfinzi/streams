"""
Utils for working with RSMs
"""

import numpy as np

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