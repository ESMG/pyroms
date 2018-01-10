import numpy as np
try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF
import pyroms


def remap(src_array, remap_file, src_grad1=None, src_grad2=None, \
             src_grad3=None, spval=1e37, verbose=False):
    '''
    remap based on addresses and weights computed in a setup phase
    '''

    # get info from remap_file
    data = netCDF.Dataset(remap_file, 'r')
    title = data.title
    map_method = data.map_method
    normalization = data.normalization
    src_grid_name = data.source_grid
    dst_grid_name = data.dest_grid
    src_grid_size = len(data.dimensions['src_grid_size'])
    dst_grid_size = len(data.dimensions['dst_grid_size'])
    num_links = len(data.dimensions['num_links'])
    src_grid_dims = data.variables['src_grid_dims'][:]
    dst_grid_dims = data.variables['dst_grid_dims'][:]

    # get weights and addresses from remap_file
    map_wts = data.variables['remap_matrix'][:]
    dst_add = data.variables['dst_address'][:]
    src_add = data.variables['src_address'][:]

    # get destination mask
    dst_mask = data.variables['dst_grid_imask'][:]

    # remap from src grid to dst grid
    if src_grad1 is not None:
        iorder = 2
    else:
        iorder = 1

    if verbose is True:
        print('Reading remapping: ', title)
        print('From file: ', remap_file)
        print(' ')
        print('Remapping between:')
        print(src_grid_name)
        print('and')
        print(dst_grid_name)
        print('Remapping method: ', map_method)

    ndim = len(src_array.squeeze().shape)

    if (ndim == 2):
        tmp_dst_array = np.zeros((dst_grid_size))
        tmp_src_array = src_array.flatten()

        if iorder == 1:
            # first order remapping
            # insure that map_wts is a (num_links,4) array
            tmp_map_wts = np.zeros((num_links,4))
            tmp_map_wts[:,0] = map_wts[:,0].copy()
            map_wts = tmp_map_wts
            pyroms.remapping.scrip.remap(tmp_dst_array, map_wts, \
                                         dst_add, src_add, tmp_src_array)

        if iorder == 2:
            # second order remapping
            if map_method == 'conservative':
                # insure that map_wts is a (num_links,4) array
                tmp_map_wts = np.zeros((num_links,4))
                tmp_map_wts[:,0:2] = map_wts[:,0:2].copy()
                map_wts = tmp_map_wts
                tmp_src_grad1 = src_grad1.flatten()
                tmp_src_grad2 = src_grad2.flatten()
                pyroms.remapping.scrip.remap(tmp_dst_array, map_wts, \
                                             dst_add, src_add, tmp_src_array, \
                                             tmp_src_grad1, tmp_src_grad2)
            elif map_method == 'bicubic':
                tmp_src_grad1 = src_grad1.flatten()
                tmp_src_grad2 = src_grad2.flatten()
                tmp_src_grad3 = src_grad3.flatten()
                pyroms.remapping.scrip.remap(tmp_dst_array, map_wts, \
                                             dst_add, src_add, tmp_src_array, \
                                             tmp_src_grad1, tmp_src_grad2, \
                                             tmp_src_grad3)
            else:
                raise ValueError('Unknown method')

        # mask dst_array
        idx = np.where(dst_mask == 0)
        tmp_dst_array[idx] = spval
        tmp_dst_array = np.ma.masked_values(tmp_dst_array, spval)

        # reshape
        dst_array = np.reshape(tmp_dst_array, (dst_grid_dims[1], \
                               dst_grid_dims[0]))

    elif (ndim == 3):

        nlev = src_array.shape[0]
        dst_array = np.zeros((nlev, dst_grid_dims[1], dst_grid_dims[0]))

        # loop over vertical level
        for k in range(nlev):

            tmp_src_array = src_array[k,:,:].flatten()
            tmp_dst_array = np.zeros((dst_grid_size))

            if iorder == 1:
                # first order remapping
                # insure that map_wts is a (num_links,4) array
                tmp_map_wts = np.zeros((num_links,4))
                tmp_map_wts[:,0] = map_wts[:,0].copy()
                map_wts = tmp_map_wts
                pyroms.remapping.scrip.remap(tmp_dst_array, map_wts, \
                                             dst_add, src_add, tmp_src_array)

            if iorder == 2:
                # second order remapping
                if map_method == 'conservative':
                    tmp_src_grad1 = src_grad1.flatten()
                    tmp_src_grad2 = src_grad2.flatten()
                    pyroms.remapping.scrip.remap(tmp_dst_array, map_wts, \
                                                 dst_add, src_add, tmp_src_array, \
                                                 tmp_src_grad1, tmp_src_grad2)
                elif map_method == 'bicubic':
                    tmp_src_grad1 = src_grad1.flatten()
                    tmp_src_grad2 = src_grad2.flatten()
                    tmp_src_grad3 = src_grad3.flatten()
                    pyroms.remapping.scrip.remap(tmp_dst_array, map_wts, \
                                                 dst_add, src_add, tmp_src_array, \
                                                 tmp_src_grad1, tmp_src_grad2, \
                                                 tmp_src_grad3)
                else:
                    raise ValueError('Unknown method')


            # mask dst_array
            idx = np.where(dst_mask == 0)
            tmp_dst_array[idx] = spval
            tmp_dst_array = np.ma.masked_values(tmp_dst_array, spval)

            # reshape
            dst_array[k,:,:] = np.reshape(tmp_dst_array, (dst_grid_dims[1], \
                                          dst_grid_dims[0]))

    else:
        raise ValueError('src_array must have two or three dimensions')


    # close data file
    data.close()

    return dst_array
