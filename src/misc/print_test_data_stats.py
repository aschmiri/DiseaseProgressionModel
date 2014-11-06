#! /usr/bin/env python2.7
from common import log as log
from vgam.datahandler import DataHandler
from common import adni_tools as adni
import os.path


def main():

    data_handler = DataHandler.get_data_handler()
    measurements = data_handler.get_measurements_as_dict(visits=['bl', 'm12', 'm24'],
                                                         biomarkers=adni.biomarker_names,
                                                         select_test_set=True,
                                                         select_complete=True)

    # Setup DB
    import rpy2.robjects as robjects
    robjects.r['load'](os.path.join(adni.merge_folder, 'data/adnimerge.rdata'))
    adnimerge = robjects.r['adnimerge']

    # Print each row
    dx_dict = {0.0: 'CN', 0.25: 'EMCI', 0.5: 'MCI', 0.75: 'LMCI', 1.0: 'AD'}
    males = {'CN': 0, 'EMCI': 0, 'LMCI': 0, 'AD': 0}
    females = {'CN': 0, 'EMCI': 0, 'LMCI': 0, 'AD': 0}
    rids = adnimerge[adnimerge.colnames.index('RID')]
    for rid in measurements:
        for row in range(adnimerge.nrow):
            if rids[row] == rid:
                gender = adnimerge[9][row]
                dx = dx_dict[measurements[rid]['bl']['DX.scan']]
                break

        if gender == 1:
            males[dx] += 1
        else:
            females[dx] += 1

    print log.RESULT, 'Males:   {0}'.format(males)
    print log.RESULT, 'Females: {0}'.format(females)

if __name__ == '__main__':
    main()
