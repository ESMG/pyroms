import datetime as dt
import netCDF4 as nc
import subprocess
#import pdb; pdb.set_trace()

class tagfile():

        def __init__(self,filename):
                self.filein = filename
                return None

        def __call__(self):
                # read time in netcdf
                time, timeunits = self.read_time()
                # create the date tag
                tag = self.create_tag(time, timeunits)
                # define a new filename
                self.create_new_filename(tag)
                # rename file
                self.rename_file()
                return None

        def read_time(self):
                ''' read ocean_time variable and units in netcdf file'''
                fid = nc.Dataset(self.filein,'r')
                time = fid.variables['ocean_time'][:]
                timeunits = fid.variables['ocean_time'].units
                fid.close()
                if len(time) > 1:
#			print 'error : multiple values in time array' ; exit()
                        ntim = len(time)
                        time = time[ntim-1]
                else:
                        time = time[0]
                return time, timeunits

        def create_tag(self,time, timeunits):
                ''' create a datetime object from reference date and ocean_time'''
                # ugly part to get reference date from units string
                units_wrk  = timeunits.replace(':',' ').replace('-',' ').split()
                delta_type = units_wrk[0]
                year_ref   = int(units_wrk[2])
                month_ref  = int(units_wrk[3])
                day_ref    = int(units_wrk[4])
                hour_ref   = int(units_wrk[5])
                min_ref    = int(units_wrk[6])
                sec_ref    = int(units_wrk[7])
                # create datetime object for reference date
                dateref = dt.datetime(year_ref,month_ref,day_ref,hour_ref,min_ref,sec_ref)
                # create a datetime object for current time
                if delta_type == 'seconds':
                        tag = dateref + dt.timedelta(seconds=time)
                return tag

        def create_new_filename(self,tag):
                ''' based on tag, generate a new filename '''
                # get rid of full path (if any)
                filein = self.filein.replace('/',' ').split()[-1]
                # get the pieces we want to keep in filename
                filein_wrk = filein.replace('_',' ').split()
                runname = filein_wrk[0]
                filetype = filein_wrk[1]
                # write our new filename
                self.fileout = runname + '_' + filetype + '_' + tag.isoformat() + '.nc'
                self.fileout = self.fileout.replace(':00:00','')
                return None


        def rename_file(self):
                ''' call unix command mv '''
                # remove filein from full path
                wrk = self.filein.replace('/',' ').split()[0:-1]
                # re-create path
                path = '.'
                for part in wrk:
                        path = path + '/' + part
                # rename file
                subprocess.call('mv ' + self.filein + ' ' + path + '/' + self.fileout, shell=True)
                return None


#----------------------------------------------------------------------------------------------
# example

lis = subprocess.check_output(["ls *_?????.nc"], shell=True)
lis = lis.split()

for file in lis:
    file = str(file).replace("b'","").replace("'","")
    print(file)
    mytag = tagfile(file)
    mytag()
