def shift_SODA_data(data):

    sdata = data.copy()

    if len(sdata.shape) == 1:
        # 1d array [lon]
        sdata[:360] = data[360:]
        sdata[360:] = data[:360]    

    if len(sdata.shape) == 2:
        # 2d array [lat,lon]
        sdata[:,:360] = data[:,360:]
        sdata[:,360:] = data[:,:360]

    if len(sdata.shape) == 3:
        # 3d array [depth,lat,lon]
        sdata[:,:,:360] = data[:,:,360:]
        sdata[:,:,360:] = data[:,:,:360]
    
    if len(sdata.shape) == 4:
        # 4d array [time,depth,lat,lon]
        sdata[:,:,:,:360] = data[:,:,:,360:]
        sdata[:,:,:,360:] = data[:,:,:,:360]

    return sdata 
