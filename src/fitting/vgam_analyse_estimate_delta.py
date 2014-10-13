#! /usr/bin/env python2.7
import argparse
import numpy as np
from common import log as log
from vgam.datahandler import DataHandler
import fitting.vgam_evaluation as ve


def main():
    parser = argparse.ArgumentParser()
    parser = DataHandler.add_arguments(parser)
    parser.add_argument('visits', nargs='+', type=str)
    parser.add_argument('--consistent_data', action='store_true', help='us only subjects with bl, m12 and m24 visits')
    parser.add_argument('--estimate_dprs', action='store_true', help='estimate dpis and dprs')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpi / dpr estimations')
    parser.add_argument('--recompute_predictions', action='store_true', help='recompute the biomarker predictions')
    parser.add_argument('--exclude_cn', action='store_true', help='exclude healthy subjects from analysis')

    estimates = {}
    methods = ['cog', 'long', 'mbl', 'img', 'all']
    for method in methods:
        estimates.update({method: {}})
        for visits in [['bl'], ['m12'], ['m24'], ['bl', 'm12'], ['m12', 'm24']]:
            visits_string = '_'.join(visits)
            estimates[method].update({visits_string: {}})
            args = parser.parse_args('{0} -m {1} --consistent_data'.format(' '.join(visits), method).split())
            _, diagnoses, dpis, _, _, _ = ve.get_progression_estimates(args)
            diagnoses = np.array(diagnoses)
            dpis = np.array(dpis)
            estimates[method][visits_string].update({'CN': np.mean(dpis[np.where(diagnoses == 0.0)])})
            estimates[method][visits_string].update({'EMCI': np.mean(dpis[np.where(diagnoses == 0.25)])})
            estimates[method][visits_string].update({'LMCI': np.mean(dpis[np.where(diagnoses == 0.75)])})
            estimates[method][visits_string].update({'AD': np.mean(dpis[np.where(diagnoses == 1.0)])})

    for method in methods:
        print log.INFO, 'Results for {0}'.format(method)
        for diagnosis in ['CN', 'EMCI', 'LMCI', 'AD']:
            # print log.RESULT, '{0: <4}:   {1:.2f} {2:.2f} {3:.2f} | {4:.2f} {5:.2f}'.format(
            #                  diagnosis,
            #                  estimates[method]['bl'][diagnosis],
            #                  estimates[method]['m12'][diagnosis],
            #                  estimates[method]['m24'][diagnosis],
            #                  estimates[method]['bl_m12'][diagnosis],
            #                  estimates[method]['m12_m24'][diagnosis])
            print log.RESULT, '{0: <4}:   {1:.2f} {2:.2f} | {3:.2f}  '.format(
                              diagnosis,
                              estimates[method]['m12'][diagnosis] - estimates[method]['bl'][diagnosis],
                              estimates[method]['m24'][diagnosis] - estimates[method]['m12'][diagnosis],
                              estimates[method]['m12_m24'][diagnosis] - estimates[method]['bl_m12'][diagnosis])

if __name__ == '__main__':
    main()
