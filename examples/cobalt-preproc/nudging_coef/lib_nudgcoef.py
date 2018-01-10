import numpy as npy
import netCDF4 as nc
import pyroms

class nudgcoef():
    ''' A class to write the Nudging coeficient file for ROMS '''

    def __init__(self,roms_grid):
        ''' init an object of the class with the pyroms grid ID '''
        self.grd = pyroms.grid.get_ROMS_grid(roms_grid)
        return None

    def __call__(self,east_dict,west_dict,north_dict,south_dict,tracer_timescales,foutname='./nudging_coef.nc'):
        ''' call with following dictionaries :
        4 boundaries dict + tracer timescales
        for example :
        east_dict  = {'nudge':True,'factor': 1,'width':50,'transition':'linear'}
        west_dict  = {'nudge':True,'factor': 1,'width':50,'transition':'linear'}
        north_dict = {'nudge':True,'factor': 1,'width':50,'transition':'linear'}
        south_dict = {'nudge':True,'factor': 1,'width':50,'transition':'linear'}
        tracer_timescales = {'M2':30,'M3':30,'temp':30,'salt':30,'tracer':30}
        tips:
        * nudge = True if open boundary, False otherwise
        * factor allows to have different timescales at each boundary
        * width is in grid points
        * transition shapes how timescale varies spatially
        * tracer timescales are in days
        '''
        self.east_dict = east_dict
        self.west_dict = west_dict
        self.north_dict = north_dict
        self.south_dict = south_dict
        self.tra_ts = tracer_timescales
        self.foutname = foutname
        # create 2d coef
        self.nud2 = self._create_nudgcoef_2d()
        # create 3d coef
        self.nud3 = self._create_nudgcoef_3d()
        # write to netcdf
        self._write_nc_file()
        return None

    def _create_nudgcoef_3d(self):
        ''' expand 2d coef along the vertical '''
        # RD: later we could imagine multiplying by
        # a vertical profile if needed
        ny, nx = self.grd.hgrid.mask_rho.shape
        nz = self.grd.vgrid.N
        nudgcoef = npy.zeros((nz,ny,nx))
        for kz in npy.arange(nz):
            nudgcoef[kz,:,:] = self.nud2[:,:]
        return nudgcoef

    def _create_nudgcoef_2d(self):
        ''' create the 2d nudging coef from dictionaries '''
        ny, nx = self.grd.hgrid.mask_rho.shape
        nudgcoef_west  = npy.zeros((ny,nx))
        nudgcoef_east  = npy.zeros((ny,nx))
        nudgcoef_north = npy.zeros((ny,nx))
        nudgcoef_south = npy.zeros((ny,nx))
        nudgcoef       = npy.zeros((ny,nx))
        mask           = self.grd.hgrid.mask_rho
        # west boundary
        if self.west_dict['nudge'] is True:
            fc = self.west_dict['factor']
            wd = self.west_dict['width']
            tr = self.west_dict['transition']
            if tr == 'linear':
                for ji in npy.arange(0,wd):
                    nudgcoef_west[:,ji] = fc * (wd-ji) / float(wd)
            elif tr == 'linear_nocoast':
                for ji in npy.arange(0,wd):
                    nudgcoef_west[:,ji] = mask[:,0] * fc * (wd-ji) / float(wd)
            else:
                print('transition not coded') ; pass
        # east boundary
        if self.east_dict['nudge'] is True:
            fc = self.east_dict['factor']
            wd = self.east_dict['width']
            tr = self.east_dict['transition']
            if tr == 'linear':
                for ji in npy.arange(nx-wd,nx):
                    nudgcoef_east[:,ji] = fc * (wd-nx+ji) / float(wd)
            elif tr == 'linear_nocoast':
                for ji in npy.arange(nx-wd,nx):
                    nudgcoef_east[:,ji] = mask[:,-1] * fc * (wd-nx+ji) / float(wd)
            else:
                print('transition not coded') ; pass
        # south boundary
        if self.south_dict['nudge'] is True:
            fc = self.south_dict['factor']
            wd = self.south_dict['width']
            tr = self.south_dict['transition']
            if tr == 'linear':
                for jj in npy.arange(0,wd):
                    nudgcoef_south[jj,:] = fc * (wd-jj) / float(wd)
            if tr == 'linear_nocoast':
                for jj in npy.arange(0,wd):
                    nudgcoef_south[jj,:] = mask[0,:] * fc * (wd-jj) / float(wd)
            else:
                print('transition not coded') ; pass
        # north boundary
        if self.north_dict['nudge'] is True:
            fc = self.north_dict['factor']
            wd = self.north_dict['width']
            tr = self.north_dict['transition']
            if tr == 'linear':
                for jj in npy.arange(ny-wd,ny):
                    nudgcoef_south[jj,:] = fc * (wd-ny+jj) / float(wd)
            if tr == 'linear_nocoast':
                for jj in npy.arange(ny-wd,ny):
                    nudgcoef_south[jj,:] = mask[-1,:] * fc * (wd-ny+jj) / float(wd)
            else:
                print('transition not coded') ; pass


        # create the total coefficient by combining all 4 fields
        # the max functions is useful to make nice corners when
        # individual field overlap
        # maybe not the most efficient but short and readable
        for jj in npy.arange(ny):
            for ji in npy.arange(nx):
                nudgcoef[jj,ji] = max(nudgcoef_west[jj,ji], \
                nudgcoef_east[jj,ji],nudgcoef_north[jj,ji],nudgcoef_south[jj,ji])
        return nudgcoef

    def _write_nc_file(self):
        ''' writing to netcdf and multiplying by inverse timescales '''
        ncfile = self.foutname
        fid = nc.Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')

            # dimensions
        fid.createDimension('xi_rho', npy.size(self.grd.hgrid.mask_rho,1))
        fid.createDimension('eta_rho', npy.size(self.grd.hgrid.mask_rho,0))
        fid.createDimension('s_rho', self.grd.vgrid.N)
        fid.createDimension('s_w', self.grd.vgrid.Np)
            fid.description = 'Nudging coefficients for grid' + self.grd.name
        # vertical coordinate
        fid.createVariable('s_rho', 'f8', ('s_rho'))
        fid.variables['s_rho'].long_name = 'S-coordinate at RHO-points'
        fid.variables['s_rho'].valid_min = '-1'
        fid.variables['s_rho'].valid_max = '0'
        fid.variables['s_rho'].field = 's_rho,scalar'
        fid.variables['s_rho'][:] = self.grd.vgrid.s_rho

            # variables
        O_M2_NudgeCoef     = fid.createVariable('M2_NudgeCoef',     'f8', ('eta_rho','xi_rho',))
        O_M3_NudgeCoef     = fid.createVariable('M3_NudgeCoef',     'f8', ('s_rho','eta_rho','xi_rho',))
        O_temp_NudgeCoef   = fid.createVariable('temp_NudgeCoef',   'f8', ('s_rho','eta_rho','xi_rho',))
        O_salt_NudgeCoef   = fid.createVariable('salt_NudgeCoef',   'f8', ('s_rho','eta_rho','xi_rho',))
        O_tracer_NudgeCoef = fid.createVariable('tracer_NudgeCoef', 'f8', ('s_rho','eta_rho','xi_rho',))
        # data
        O_M2_NudgeCoef[:,:]       = (1./self.tra_ts['M2'])     * self.nud2
        O_M3_NudgeCoef[:,:,:]     = (1./self.tra_ts['M3'])     * self.nud3
        O_temp_NudgeCoef[:,:,:]   = (1./self.tra_ts['temp'])   * self.nud3
        O_salt_NudgeCoef[:,:,:]   = (1./self.tra_ts['salt'])   * self.nud3
        O_tracer_NudgeCoef[:,:,:] = (1./self.tra_ts['tracer']) * self.nud3
        # attributes
        O_M2_NudgeCoef.long_name = '2D momentum inverse nudging coefficients'
        O_M2_NudgeCoef.units = 'days-1'
        O_M2_NudgeCoef.coordinates = 'xi_rho eta_rho'

        O_M3_NudgeCoef.long_name = '3D momentum inverse nudging coefficients'
        O_M3_NudgeCoef.units = 'days-1'
        O_M3_NudgeCoef.coordinates = 'xi_rho eta_rho s_rho'

        O_temp_NudgeCoef.long_name = 'temp inverse nudging coefficients'
        O_temp_NudgeCoef.units = 'days-1'
        O_temp_NudgeCoef.coordinates = 'xi_rho eta_rho s_rho'

        O_salt_NudgeCoef.long_name = 'salt inverse nudging coefficients'
        O_salt_NudgeCoef.units = 'days-1'
        O_salt_NudgeCoef.coordinates = 'xi_rho eta_rho s_rho'

        O_tracer_NudgeCoef.long_name = 'generic tracer inverse nudging coefficients'
        O_tracer_NudgeCoef.units = 'days-1'
        O_tracer_NudgeCoef.coordinates = 'xi_rho eta_rho s_rho'
        # close
        fid.close()
        return None


#----------------------------------------------------------------------------
# example :

ccs1 = nudgcoef('CCS')

east  = {'nudge':False,'factor': 1,'width':10,'transition':'linear_nocoast'}
west  = {'nudge':True,'factor': 1,'width':10,'transition':'linear_nocoast'}
north = {'nudge':True,'factor': 1,'width':10,'transition':'linear_nocoast'}
south = {'nudge':True,'factor': 1,'width':10,'transition':'linear_nocoast'}

tracer_timescales = {'M2':30,'M3':30,'temp':30,'salt':30,'tracer':30}

#ccs1(east,west,north,south,tracer_timescales)

# strong restoring
east  = {'nudge':False,'factor': 1,'width':70,'transition':'linear_nocoast'}
west  = {'nudge':True,'factor': 1,'width':70,'transition':'linear_nocoast'}
north = {'nudge':True,'factor': 1,'width':70,'transition':'linear_nocoast'}
south = {'nudge':True,'factor': 1,'width':70,'transition':'linear_nocoast'}

tracer_timescales = {'M2':10,'M3':10,'temp':10,'salt':10,'tracer':10}

ccs1(east,west,north,south,tracer_timescales,foutname='./CCS_nudging_coef_large.nc')
