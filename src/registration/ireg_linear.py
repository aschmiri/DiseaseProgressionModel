#! /usr/bin/env python

EXEC_REG = 'ireg'
EXEC_TRANS = 'transformation'


def run(source, template, mask, dof_file, rreg_param_file, areg_param_file, warped_img=None, verbose=True):
    """
    Performs a linear registration using the IRTK ireg program.
    """
    import os.path
    from subprocess import call

    dofRigid = dof_file.replace('.dof.gz', '_rigid.dof.gz')
    if os.path.exists(dof_file):
        if verbose:
            print 'File', dof_file, 'already exists'
    else:
        if os.path.exists(dofRigid):
            if verbose:
                print 'File', dof_file, 'already exists'
        else:
            #
            # Run rigid registration
            #
            if verbose:
                print '--------------------'
                print 'Starting: ireg'
                print 'Template:', template
                print 'Source:  ', source
                print 'Mask:    ', mask
                print 'DOF out: ', dofRigid
                print 'Param:   ', rreg_param_file

            if mask in [None, 'None', 'none']:
                call([EXEC_REG, template, source, '-dofout', dofRigid, '-parin', rreg_param_file])
            else:
                call([EXEC_REG, template, source, '-dofout', dofRigid, '-parin', rreg_param_file, '-mask', mask])

        #
        # Run affine registration
        #
        if verbose:
            print '--------------------'
            print 'Starting: ireg'
            print 'Template:', template
            print 'Source:  ', source
            print 'Mask:    ', mask
            print 'DOF in:  ', dofRigid
            print 'DOF out: ', dof_file
            print 'Param:   ', areg_param_file

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
                print 'File', warped_img, 'already exists'
        else:
            if verbose:
                print '--------------------'
                print 'Starting: transformation'
                print 'Source:  ', source
                print 'Template:', template
                print 'DOF:     ', dof_file
                print 'Deformed:', warped_img

            call([EXEC_TRANS, source, warped_img, '-target', template, '-dofin', dof_file, '-cspline', '-matchInputType', '-Sp', '0'])
