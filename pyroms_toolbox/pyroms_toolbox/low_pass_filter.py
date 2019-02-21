import numpy as np

def low_pass_filter(data, window_size=3):
    """
    A simple low pass filter using a square window moving average

    Parameters
    ----------
    data : array_like
        input array. Can be multi-dimensional but time need to be the first dimension.
    window_size : array_like, optional
        A scalar giving the size of the filter window. window_size must be odd.
    
    Returns
    -------
    low_passed_data : ndarray
        An array the same size as input containing the low-passed result.
    """

    low_passed_data = np.zeros(data.shape)

    if np.abs(window_size) % 2 != 1:
        raise ValueError('window_size must be odd.')

    win2 = (window_size - 1) / 2

    if (len(data.squeeze().shape) == 1):
        for n in range(np.size(data,axis=0)):
            idx1 = int(np.maximum(0, n-win2))
            idx2 = int(np.minimum(np.size(data,axis=0), n+win2+1))
            low_passed_data[n] = data[idx1:idx2].sum(axis=0) / float(idx2-idx1)

    else:
        for n in range(np.size(data,axis=0)):
            idx1 = int(np.maximum(0, n-win2))
            idx2 = int(np.minimum(np.size(data,axis=0), n+win2+1))
            low_passed_data[n,:] = data[idx1:idx2,:].sum(axis=0) / float(idx2-idx1)

    return low_passed_data
