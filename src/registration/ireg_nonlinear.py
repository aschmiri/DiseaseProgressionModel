#! /usr/bin/env python

def run( source, template, dofIn, dofOut, paramFile, wapedImg ):
    """
    Performs a nonlinear registration using the IRTK ireg program (currently called reg3).
    """
    import os.path
    from subprocess import call
    
    iregExec = 'ireg'
    transExec = 'transformation'
    
    #
    # Run nonlinear registration
    #
    if os.path.exists( dofOut ):
        print 'File', dofOut, 'already exists'
    else:
        print '--------------------'
        print 'Starting: ireg'
        print 'Template: ', template
        print 'Source:   ', source
        print 'DOF in:   ', dofIn
        print 'DOF out:  ', dofOut
        print 'Param:    ', paramFile
        
        if dofIn in [None, 'None', 'none']:
            call([ iregExec, template, source, '-dofout', dofOut, '-parin', paramFile, '-v' ])
        else:
            call([ iregExec, template, source, '-dofin', dofIn, '-dofout', dofOut, '-parin', paramFile, '-v' ])

    
    #
    # Run transformation
    #
    if not wapedImg in [None, 'None', 'none']:
        if os.path.exists( wapedImg ):
            print 'File', wapedImg, 'already exists'
        else:
            print '--------------------'
            print 'Starting transformation'
            print 'Source:  ', source
            print 'Template:', template
            print 'DOF:     ', dofOut
            print 'Deformed:', wapedImg
            
            call([ transExec, source, wapedImg, '-target', template, '-dofin', dofOut, '-cspline', '-matchInputType', '-Sp', '0' ])
