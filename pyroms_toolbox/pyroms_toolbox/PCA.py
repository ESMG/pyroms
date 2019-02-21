import numpy as np
from matplotlib.mlab import psd
from scipy.signal import get_window 



class PCA(object):

    def __init__(self, dataset, algorithm='svd'):
        """
        This class carries out Principal Component Analysis on dataset.

        Arguments:
            'dataset' -- NumPy array containing the field to be decomposed. 
                         First dimension must be time.

        The optional algorithm parameter can be either 'svd' to perform PCA with 
        the singular value decomposition, or 'eig' to use a symmetric eigenvalue
        decomposition. 
        
        Notes:
          PCA does not center or scale dataset; you usually want to first
            dataset = center(dataset) or
            dataset = standardize(dataset)
          with the functions provided.
        """

        self.dataset = dataset

        self.packedata, self.mask = self.pack(self.dataset)

        if algorithm == 'svd':
            self.packedEOFs, self.ECs, self.L = self._pca_svd(self.packedata)
        elif algorithm == 'eig':
            self.packedEOFs, self.ECs, self.L = self._pca_eig(self.packedata)
        else:
            raise RuntimeError('Algorithm %s not known.'%algorithm)

        self.EOFs = self.unpack(self.packedEOFs, self.mask)

    def GetEOFs(self):
        """
        Returns the Empirical Orthogonal Functions EOFs.
        """
        return self.EOFs

    def GetECs(self):
        """
        Returns the Temporal Expansion Coeficients ECs.
        """
        return self.ECs

    def GetL(self):
        """
        Returns Covariances matrix L.
        """
        return self.L

    def GetPSD(self, blocks_length=0, fs=1, window='boxcar', overlap=0):
        """
        Returns the power spectral density estimated by Welch's average 
        periodogram method.
        ECs are divided into blocks of lenght "blocks_length". If "blocks_length = 0",
        only one block spanning the whole time serie is used. Each block is windowed 
        by the function "window". "overlap" gives the length of the overlap between 
        blocks. The PSD of each segment are then averaged.
        """
        nECs, nt = self.ECs.shape
        if blocks_length == 0:
            blocks_length = nt
        else:
            blocks_length = blocks_length
        window_vector = get_window(window, blocks_length)
        PSD = np.zeros((nECs, blocks_length/2+1))
        freq = np.zeros((nECs, blocks_length/2+1))
        for n in range(nECs):
            P, f = psd(self.ECs[n], NFFT=blocks_length, Fs=fs, window=window_vector, noverlap=overlap)
            PSD[n] = np.squeeze(P)
            freq[n] = f
        return PSD, freq
        
    def pack(self, dataset):
        """
        """
        nt = dataset.shape[0]
        mask = np.ma.getmaskarray(dataset[0])
        npt = np.sum(~mask)
        packedata = np.zeros((nt,npt))

        for t in range(nt):
            packedata[t] = np.ma.compressed(dataset[t])

        return packedata, mask


    def unpack(self, packedeofs, mask):
        """
        """
        neof = packedeofs.shape[1]
        shape = mask.shape

        if len(mask.shape) == 2:
            ny, nx = mask.shape
            dims = nx * ny
        elif len(mask.shape) == 3:
            nz, ny, nx = mask.shape
            dims = nx * ny * nz

        eofs = np.ma.masked_all((neof,dims))

        for n in range(neof):
            eofs[n,~mask.flatten()] = packedeofs[:,n]

        if len(mask.shape) == 2:
            eofs = eofs.reshape(neof,ny,nx)
        elif len(mask.shape) == 3:
            eofs = eofs.reshape(neof,nz,ny,nx)

        return eofs


    def _pca_svd(self, packedata):
        """
        Calculates EOF decomposition of a field by means 
        of the SVD decomposition.

        Returns EOFs, ECs, eigenvalues, eigenvectors


        Use Singular Value Decomposition (SVD) to break up
        M(mxn) into 3 matrices:

            M = U * S * V.T

        where U and V are orthonormal and D is diagonal.  Then,
            EOFs = V
            ECs = U * S
            eigenvalues = ECs.T * ECs / (n-1) = S**2 / (n-1)
        """
        u, s, vt = np.linalg.svd(packedata, full_matrices = False)
        packedpcs = np.transpose(vt)
        ecs = u * s
        ecs = ecs.T
        eigenval = s**2 / (len(s)-1)
        return packedpcs, ecs, eigenval


    def _pca_eig(self, packedata):
        """
        Calculates EOF decomposition of a field by means
        of the symmetric eigenvalue decomposition method
        """
        vals, vecs = self._sym_eigh(packedata)
        packedpcs = vecs
        ecs = np.dot(packedata, vecs)
        ecs = ecs.T
        eigenval = vals / (len(vals)-1)
        return packedpcs, ecs, eigenval


    def _sym_eigh(self, a):
        """
        Return the eigenvectors and eigenvalues of the symmetric matrix a'a. If
        a has more columns than rows, then that matrix will be rank-deficient,
        and the non-zero eigenvalues and eigenvectors can be more easily extracted
        from the matrix aa'.
        From the properties of the SVD:
          if a of shape (m,n) has SVD u*s*v', then:
              a'a = v*s's*v'
              aa' = u*ss'*u'
          let s_hat, an array of shape (m,n), be such that s * s_hat = I(m,m) 
          and s_hat * s = I(n,n). (Note that s_hat is just the elementwise 
          reciprocal of s, as s is zero except on the main diagonal.)

          Thus, we can solve for u or v in terms of the other:
              v = a'*u*s_hat'
              u = a*v*s_hat
        """
        m, n = a.shape
        if m >= n:
            # just return the eigenvalues and eigenvectors of a'a
            vals, vecs = self._eigh(np.dot(np.transpose(a), a))
            vecs = np.where(vecs < 0, 0, vecs)
            return vals, vecs
        else:
            # figure out the eigenvalues and vectors based on aa', which is smaller
            w, v = self._eigh(np.dot(a, a.transpose()))
            # in case due to numerical instabilities we have w < 0 anywhere, 
            # peg them to zero
            vals = np.where(w < 0, 0, w)
            # now get the inverse square root of the diagonal, which will form the
            # main diagonal of s_hat
            err = np.seterr(divide='ignore', invalid='ignore')
            s_hat = 1/np.sqrt(w)
            np.seterr(**err)
            s_hat[~np.isfinite(s_hat)] = 0
            # s_hat is a list of length m, a'u is (n,m), so we can just use
            # numpy's broadcasting instead of matrix multiplication, and only create
            # the upper mxm block of a'u, since that's all we'll use anyway...
            vecs = np.dot(np.transpose(a), v[:,:m]) * s_hat
            return vals, vecs

    def _eigh(self, a):
        vals, vecs = np.linalg.eigh(a)
        order = np.flipud(vals.argsort())
        return vals[order], vecs[:,order]





def center(dataset):
    """
    Returns a centered version (mean along _first_ axis removed) of an array
    """
    return dataset - dataset.mean(axis=0)


def standardize(dataset):
    """
    Returns a standardized (centered and unit variance) of an array
    """
    residual = center(dataset)
    return residual / residual.std(axis=0)

