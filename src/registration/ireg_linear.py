#! /usr/bin/env python
import os.path
from subprocess import call
from src.common import adni_tools as adni

EXEC_REG = 'ireg'
EXEC_TRANS = 'transformation'


def run(source, template, mask, dof_file, rreg_param_file, areg_param_file, warped_img=None, verbose=True):
    '''
    Perform a linear registration using the IRTK ireg program.
    '''
    dofRigid = dof_file.replace('.dof.gz', '_rigid.dof.gz')
    if os.path.exists(dof_file):
        if verbose:
            print adni.SKIP, 'File', dof_file, 'already exists'
    else:
        if os.path.exists(dofRigid):
            if verbose:
                print adni.SKIP, 'File', dof_file, 'already exists'
        else:
            #
            # Run rigid registration
            #
            if verbose:
                print adni.INFO, '--------------------'
                print adni.INFO, 'Starting: ireg'
                print adni.INFO, 'Template:', template
                print adni.INFO, 'Source:  ', source
                print adni.INFO, 'Mask:    ', mask
                print adni.INFO, 'DOF out: ', dofRigid
                print adni.INFO, 'Param:   ', rreg_param_file

            if mask in [None, 'None', 'none']:
                call([EXEC_REG, template, source, '-dofout', dofRigid, '-parin', rreg_param_file])
            else:
                call([EXEC_REG, template, source, '-dofout', dofRigid, '-parin', rreg_param_file, '-mask', mask])

        #
        # Run affine registration
        #
        if verbose:
            print adni.INFO, '--------------------'
            print adni.INFO, 'Starting: ireg'
            print adni.INFO, 'Template:', template
            print adni.INFO, 'Source:  ', source
            print adni.INFO, 'Mask:    ', mask
            print adni.INFO, 'DOF in:  ', dofRigid
            print adni.INFO, 'DOF out: ', dof_file
            print adni.INFO, 'Param:   ', areg_param_file

        if mask in [None, 'None', 'none']:
            call([EXEC_REG, template, source, '-dofin', dofRigid, '-dofout', dof_file, '-parin', areg_param_file])
        else:
            call([EXEC_REG, template, source, '-dofin', dofRigid, '-dofout', dof_file, '-parin', areg_param_file, '-mask', mask])

        #
        # Delete temporary file
        #
        if os.path.exists(dofRigid):
            call(['rm', dofRigid])

    #
    # Run transformation
    #
    if warped_img not in [None, 'None', 'none']:
        if os.path.exists(warped_img):
            if verbose:
                print adni.SKIP, 'File', warped_img, 'already exists'
        else:
            if verbose:
                print adni.INFO, '--------------------'
                print adni.INFO, 'Starting: transformation'
                print adni.INFO, 'Source:  ', source
                print adni.INFO, 'Template:', template
                print adni.INFO, 'DOF:     ', dof_file
                print adni.INFO, 'Deformed:', warped_img

            call([EXEC_TRANS, source, warped_img, '-target', template, '-dofin', dof_file, '-cspline', '-matchInputType', '-Sp', '0'])
