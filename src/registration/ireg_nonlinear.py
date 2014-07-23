#! /usr/bin/env python

EXEC_REG = 'ireg'
EXEC_TRANS = 'transformation'


def run(source, template, dofIn, dofOut, paramFile, wapedImg, verbose=True):
    """
    Performs a nonlinear registration using the IRTK ireg program (currently called reg3).
    """
    import os.path
    from subprocess import call

    #
    # Run nonlinear registration
    #
    if os.path.exists(dofOut):
        if verbose:
            print 'File', dofOut, 'already exists'
    else:
        if verbose:
            print '--------------------'
            print 'Starting: ireg'
            print 'Template: ', template
            print 'Source:   ', source
            print 'DOF in:   ', dofIn
            print 'DOF out:  ', dofOut
            print 'Param:    ', paramFile

        if dofIn in [None, 'None', 'none']:
            call([EXEC_REG, template, source, '-dofout', dofOut, '-parin', paramFile, '-v'])
        else:
            call([EXEC_REG, template, source, '-dofin', dofIn, '-dofout', dofOut, '-parin', paramFile, '-v'])

    #
    # Run transformation
    #
    if not wapedImg in [None, 'None', 'none']:
        if os.path.exists(wapedImg):
            if verbose:
                print 'File', wapedImg, 'already exists'
        else:
            if verbose:
                print '--------------------'
                print 'Starting transformation'
                print 'Source:  ', source
                print 'Template:', template
                print 'DOF:     ', dofOut
                print 'Deformed:', wapedImg

            call([EXEC_TRANS, source, wapedImg, '-target', template, '-dofin', dofOut, '-cspline', '-matchInputType', '-Sp', '0'])
