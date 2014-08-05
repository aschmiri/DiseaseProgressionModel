#! /usr/bin/env python

EXEC_REG = 'ireg'
EXEC_TRANS = 'transformation'


def run(source, template, dof_in, dof_out, param_file, waped_img=None, verbose=True):
    """
    Performs a nonlinear registration using the IRTK ireg program (currently called reg3).
    """
    import os.path
    from subprocess import call

    #
    # Run nonlinear registration
    #
    if os.path.exists(dof_out):
        if verbose:
            print 'File', dof_out, 'already exists'
    else:
        if verbose:
            print '--------------------'
            print 'Starting: ireg'
            print 'Template: ', template
            print 'Source:   ', source
            print 'DOF in:   ', dof_in
            print 'DOF out:  ', dof_out
            print 'Param:    ', param_file

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
                print 'File', waped_img, 'already exists'
        else:
            if verbose:
                print '--------------------'
                print 'Starting transformation'
                print 'Source:  ', source
                print 'Template:', template
                print 'DOF:     ', dof_out
                print 'Deformed:', waped_img

            call([EXEC_TRANS, source, waped_img, '-target', template, '-dofin', dof_out, '-cspline', '-matchInputType', '-Sp', '0'])
