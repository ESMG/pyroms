import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
import subprocess
import numpy as np
from multiprocessing import Pool

import pyroms
import pyroms_toolbox

# This script works better sending stdout to a file, then executing that.

year = sys.argv[1]

data_dir = './data/'

dst_dir='./months/'

lst = subprocess.getoutput('ls ' + data_dir + 'Pacific_*' + year + '.nc')
lst_file = lst.split()

#print 'Make monthly average files from the following file list:'
#print lst_file
#print ' '

def do_file(file):
    index = file.find('Pacific')
    outfile = dst_dir + file[index:]
#    print 'Input filename ', file
#    print 'Output filename ', outfile, index

    command = ('ncra', '-d', 'ocean_time,0,9', file, outfile+'.01')
    print('ncra -d ocean_time,0,9', file, outfile+'.01')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,10,19', file, outfile+'.02')
    print('ncra -d ocean_time,10,19', file, outfile+'.02')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,20,29', file, outfile+'.03')
    print('ncra -d ocean_time,20,29', file, outfile+'.03')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,30,39', file, outfile+'.04')
    print('ncra -d ocean_time,30,39', file, outfile+'.04')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,40,49', file, outfile+'.05')
    print('ncra -d ocean_time,40,49', file, outfile+'.05')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,50,59', file, outfile+'.06')
    print('ncra -d ocean_time,50,59', file, outfile+'.06')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,60,69', file, outfile+'.07')
    print('ncra -d ocean_time,60,69', file, outfile+'.07')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,70,79', file, outfile+'.08')
    print('ncra -d ocean_time,70,79', file, outfile+'.08')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,80,89', file, outfile+'.09')
    print('ncra -d ocean_time,80,89', file, outfile+'.09')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,90,99', file, outfile+'.10')
    print('ncra -d ocean_time,90,99', file, outfile+'.10')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,100,109', file, outfile+'.11')
    print('ncra -d ocean_time,100,109', file, outfile+'.11')
#    subprocess.check_call(command)
    command = ('ncra', '-d', 'ocean_time,110,119', file, outfile+'.12')
    print('ncra -d ocean_time,110,119', file, outfile+'.12')
#    subprocess.check_call(command)
    command = ('ncrcat', outfile+'.*', outfile)
    print('ncrcat', outfile+'.*', outfile)
#    subprocess.check_call(command)
    command = ('rm', outfile+'.*')
    print('rm', outfile+'.*')
#    subprocess.check_call(command)

processes = 1
p = Pool(processes)
results = p.map(do_file, lst_file)
