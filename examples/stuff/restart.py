import subprocess
import os.path
import netCDF4

num = input("Restart number? ")
file_path = "$ARCHIVE/Arctic2/run45/restart_" + num
if os.path.exists(file_path):
    print("directory exists already:", file_path)
    exit()
else:
    cmd = "mkdir " + file_path
    subprocess.call([cmd], shell=True)

cmd = "mv -i arctic2_sta.nc arctic2_sta_" + num + ".nc"
subprocess.call([cmd], shell=True)
fh = netCDF4.Dataset("arctic2_flt.nc", "r")
nt = len(fh.dimensions['ocean_time'])
nt = str(nt-1)
print(nt)
fh.close()
cmd = "mv -i arctic2_flt.nc arctic2_flt_" + num + ".nc"
subprocess.call([cmd], shell=True)
cmd = "ncks -d ocean_time,"+nt+","+nt+" arctic2_flt_" + num + ".nc arctic2_flt.nc"
subprocess.call([cmd], shell=True)
cmd = "/usr/bin/rcp *rst* bigdipper:" + file_path
#cmd = "/usr/bin/rcp *fil* *rst* bigdipper:" + file_path
subprocess.call([cmd], shell=True)
#vi ocean_arctic2.in
