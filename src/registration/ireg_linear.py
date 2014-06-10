#! /usr/bin/env python

def run( source, template, mask, dofFile, rregParamFile, aregParamFile, wapedImg, verbose=True ):
    """
    Performs a linear registration using the IRTK ireg program (currently called reg3).
    """
    import os.path
    from subprocess import call
    
    regExec = 'ireg'
    transExec = 'transformation'
    
    dofRigid = dofFile.replace( '.dof.gz', '_rigid.dof.gz' )
    
    if os.path.exists( dofFile ):
        if verbose:
            print 'File', dofFile, 'already exists'
    else:
        if os.path.exists( dofRigid ):
            if verbose:
                print 'File', dofFile, 'already exists'
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
                print 'Param:   ', rregParamFile
            
            if mask in [None, 'None', 'none']:
                call([ regExec, template, source, '-dofout', dofRigid, '-parin', rregParamFile ])
            else:
                call([ regExec, template, source, '-dofout', dofRigid, '-parin', rregParamFile, '-mask', mask ])
    
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
            print 'DOF out: ', dofFile
            print 'Param:   ', aregParamFile
        
        if mask in [None, 'None', 'none']:
            call([ regExec, template, source, '-dofin', dofRigid,  '-dofout', dofFile, '-parin', aregParamFile ])
        else:
            call([ regExec, template, source, '-dofin', dofRigid,  '-dofout', dofFile, '-parin', aregParamFile, '-mask', mask ])
        
        #
        # Delete temporary file
        #
        if os.path.exists( dofRigid ):
            call([ 'rm', dofRigid ])
    
    #
    # Run transformation
    #
    if not wapedImg in [None, 'None', 'none']:
        if os.path.exists( wapedImg ):
            if verbose:
                print 'File', wapedImg, 'already exists'
        else:
            if verbose:
                print '--------------------'
                print 'Starting: transformation'
                print 'Source:  ', source
                print 'Template:', template
                print 'DOF:     ', dofFile
                print 'Deformed:', wapedImg
          
            call([ transExec, source, wapedImg, '-target', template, '-dofin', dofFile, '-cspline', '-matchInputType', '-Sp', '0' ])
        