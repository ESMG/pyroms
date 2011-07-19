import os
import pyroms


def test_remap_weights(field_choice, interp_file, output_file): 
    '''
    test addresses and weights computed in a setup phase
    '''

    # write test namelist file
    f = open('test_remap_weights_in','w')

    f.write('&remap_inputs' + '\n')
    f.write('    field_choice = ' + str(field_choice) + '\n')
    f.write('    interp_file = \'' + str(interp_file) + '\'\n')
    f.write('    output_file = \'' + str(output_file) + '\'\n')
    f.write('/')

    f.close()

    # run test weights
    pyroms.remapping.scrip.test_remap_weights('test_remap_weights_in')

    # clean
    os.remove('test_remap_weights_in') 
