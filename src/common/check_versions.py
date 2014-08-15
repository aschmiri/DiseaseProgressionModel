#! /usr/bin/env python2.7
import sys

################################################################################
#
# main
#
################################################################################
if __name__ == "__main__":
    print 'Python:       ', sys.version

    try:
        import numpy
        print 'NumPy:        ', numpy.__version__
    except:
        print 'NumPy:         NOT INSTALLED'

    try:
        import scipy
        print 'SciPy:        ', scipy.__version__
    except:
        print 'NumPy:         NOT INSTALLED'

    try:
        import matplotlib
        print 'MatPlotLib:   ', matplotlib.__version__
    except:
        print 'MatPlotLib:    NOT INSTALLED'

    try:
        import sklearn
        print 'scikit-learn: ', sklearn.__version__
    except:
        print 'scikit-learn:  NOT INSTALLED'

    try:
        import skimage
        print 'scikit-image: ', skimage.__version__
    except:
        print 'scikit-image:  NOT INSTALLED'

    try:
        import rpy2
        print 'RPy2:         ', rpy2.__version__
    except:
        print 'RPy2:          NOT INSTALLED'

    try:
        import joblib
        print 'joblib:       ', joblib.__version__
    except:
        print 'joblib:        NOT INSTALLED'
