#! /usr/bin/env python
# print __doc__

import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import common.adni_tools as adni

def strip_lists( list1, list2 ):
    rlist1 = []
    rlist2 = []
    for i in range(len(list1)):
        if list1[i] != None and list2[i] != None:
            rlist1.append( list1[i] )
            rlist2.append( list2[i] )
    return rlist1, rlist2

years, dpi, cdrsb, adas11, adas13, faq, mmse, moca = adni.get_years_after_symptoms()

print 'Pearson Correlation Coefficients'
print 'DPI:     ', pearsonr( *strip_lists(years, dpi) )
print 'CDR-SoB: ', pearsonr( *strip_lists(years, cdrsb) )
print 'ADAS-11: ', pearsonr( *strip_lists(years, adas11) )
print 'ADAS-13: ', pearsonr( *strip_lists(years, adas13) )
print 'FAQ:     ', pearsonr( *strip_lists(years, faq) )
print 'MMSE:    ', pearsonr( *strip_lists(years, mmse) )
print 'MOCA:    ', pearsonr( *strip_lists(years, moca) )

f, ((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8)) = plt.subplots(nrows=2,ncols=4)
ax1.scatter( dpi,  dpi,  color='red' )
ax2.scatter( dpi,  cdrsb,  color='red' )
ax3.scatter( dpi,  adas11,  color='red' )
ax4.scatter( dpi,  adas13,  color='red' )
ax5.scatter( dpi,  faq,  color='red' )
ax6.scatter( dpi,  mmse,  color='red' )
ax7.scatter( dpi,  moca,  color='red' )

plt.show()
