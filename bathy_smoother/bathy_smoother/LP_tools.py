import numpy as np
from numpy.random import random
from time import localtime
import os
import logging

from lpsolve55 import *

# This code is adapted from the matlab code
# "LP Bathymetry" by Mathieu Dutour Sikiric
# http://drobilica.irb.hr/~mathieu/Bathymetry/index.html
# For a description of the method, see
# M. Dutour Sikiric, I. Janekovic, M. Kuzmic, A new approach to
# bathymetry smoothing in sigma-coordinate ocean models, Ocean
# Modelling 29 (2009) 128--136.

def WriteLinearProgram(FileName, iList, jList, sList, Constant, ObjectiveFct):

    # Do not remove the test of feasibility !
    # It is a preliminary check and it is actually useful.

    nbVar = ObjectiveFct.shape[0]
    nbConst = Constant.shape[0]
    print('Write linear program')
    print('nbvar = ', nbVar, ' nbConst = ', nbConst)
    print(' ')

    f = open(FileName,'w')

    # write the function to minimise
    # we want to minimize the perturbation P = sum_e(|dh(e)|
    # as the absolute value is not a linear function, we use a trick:
    # we introduce an additional variable ad(e) satisfying +/- dh(e) <= ad(e)
    # and we minimize sum_e(ad(e))
    f.write('min: ')
    for iVar in range(nbVar):
        eVal = ObjectiveFct[iVar,0]
        if (eVal != 0):
            if (eVal > 0):
                add='+'
            else:
                add=''

            string = '%s%f X%d ' %(add,eVal,iVar+1)
            f.write(string)
        else:
            string = '+0 X%d ' %(iVar+1)
            f.write(string)

    f.write(';\n')
    f.write('\n')

    #write the inequality constraintes
    tolCrit = 1e-6
    for iConst in range(nbConst):
        H = np.where(iList == iConst+1)[0]
        nbH = H.shape[0]
        if (nbH == 0):
            if (Constant[iConst,0] < -tolCrit):
                testfeasibility = 0
                return testfeasibility

        else:
            string = 'row%s: ' %str(iConst+1)
            f.write(string)
            for iH in range(nbH):
                jL = jList[H[iH],0]
                sL = sList[H[iH],0]
                string = '%.2f X%d ' %(sL, jL)
                if (sL > 0):
                    add='+'
                else:
                    add=''

                string = '%s%s' %(add, string)
                f.write(string)

            string = '<= %.8e ;\n' %Constant[iConst,0]
            f.write(string)

    f.write('\n')

    # the free command does not seem to work as advertised
    f.write('free')
    for iVar in range(nbVar):
        if (iVar+1 > 1):
            f.write(',')

        string = ' X%d' %(iVar+1)
        f.write(string)

    f.write(';\n')
    f.close()
    testfeasibility = 1


    return testfeasibility



def SolveLinearProgram(iList, jList, sList, Constant, ObjectiveFct):

    nbVar = ObjectiveFct.shape[0]
    nbConstraint = Constant.shape[0]
    print('Solving a linear program of ', nbVar, ' variables and ', nbConstraint, ' Constraints')

    while(1):
        H = localtime()
        V0 = H[3]
        V1 = H[4]
        V2 = np.floor(H[5])
        V3 = np.ceil(10 * random())
        Prefix = '/tmp/lp_%s_%s_%s_%s' %(str(V0),str(V1),str(V2),str(V3))

        FileInput = '%s_input.lp' %Prefix
        FileOutput = '%s_output.lp' %Prefix
        if (os.path.exists(FileInput) is False and os.path.exists(FileOutput) is False):
            break

        print('We failed with FileInput = ', FileInput)

    testfeasibility = WriteLinearProgram(FileInput, iList, jList, sList, Constant, ObjectiveFct)
    if (testfeasibility == 0):
        raise ValueError('Feasibility test failed. testfeasibility = 0.')

    print('Linear program written in FileOutput=', FileOutput)

    lp_handle = lpsolve('read_lp_file', FileInput)
    result = lpsolve('solve', lp_handle)
    obj, ValueVar, ValueFct, testfeasibility = lpsolve('get_solution', result)
    lpsolve('delete_lp', result)

    print('Linear program solved')


    return ValueFct, ValueVar, testfeasibility
