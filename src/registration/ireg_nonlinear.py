#! /usr/bin/env python
import os.path
from subprocess import call
from src.common import adni_tools as adni

EXEC_REG = 'ireg'
EXEC_TRANS = 'transformation'


def run(source, template, dof_in, dof_out, param_file, waped_img=None, verbose=True):
    '''
    Performs a nonlinear registration using the IRTK ireg program (currently called reg3).
    '''
    #
    # Run nonlinear registration
    #
    if os.path.exists(dof_out):
        if verbose:
            print adni.SKIP, 'File', dof_out, 'already exists'
    else:
        if verbose:
            print adni.INFO, '--------------------'
            print adni.INFO, 'Starting: ireg'
            print adni.INFO, 'Template: ', template
            print adni.INFO, 'Source:   ', source
            print adni.INFO, 'DOF in:   ', dof_in
            print adni.INFO, 'DOF out:  ', dof_out
            print adni.INFO, 'Param:    ', param_file

        if dof_in in [None, 'None', 'none']:
            call([EXEC_REG, template, source, '-dofout', dof_out, '-parin', param_file, '-v'])
        else:
            call([EXEC_REG, template, source, '-dofin', dof_in, '-dofout', dof_out, '-parin', param_file, '-v'])

    #
    # Run transformation
    #
    if waped_img not in [None, 'None', 'none']:
        if os.path.exists(waped_img):
            if verbose:
                print adni.SKIP, 'File', waped_img, 'already exists'
        else:
            if verbose:
                print adni.INFO, '--------------------'
                print adni.INFO, 'Starting transformation'
                print adni.INFO, 'Source:  ', source
                print adni.INFO, 'Template:', template
                print adni.INFO, 'DOF:     ', dof_out
                print adni.INFO, 'Deformed:', waped_img

            call([EXEC_TRANS, source, waped_img, '-target', template, '-dofin', dof_out, '-cspline', '-matchInputType', '-Sp', '0'])
