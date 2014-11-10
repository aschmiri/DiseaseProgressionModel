#! /usr/bin/env python2.7
import argparse
import numpy as np
from common import log as log
from common import evaluation_tools as et


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--consistent_data', action='store_true', help='us only subjects with bl, m12 and m24 visits')
    parser.add_argument('--estimate_dprs', action='store_true', help='estimate dpis and dprs')
    parser.add_argument('--recompute_estimates', action='store_true', help='recompute the dpi / dpr estimations')
    parser.add_argument('--recompute_predictions', action='store_true', help='recompute the biomarker predictions')
    parser.add_argument('--exclude_cn', action='store_true', help='exclude healthy subjects from analysis')
    args = parser.parse_args()

    estimates = {}
    methods = ['cog', 'vol', 'ml', 'img', 'all']
    for method in methods:
        estimates.update({method: {}})
        for visits in [['bl'], ['m12'], ['m24'], ['bl', 'm12'], ['m12', 'm24']]:
            _, diagnoses, dpis, _, _, _ = et.get_progress_estimates(visits, method=method, estimate_dprs=args.estimate_dprs)

            diagnoses = np.array(diagnoses)
            dpis = np.array(dpis)
            visits_string = '_'.join(visits)
            estimates[method].update({visits_string: {}})
            estimates[method][visits_string].update({'CN': np.mean(dpis[np.where(diagnoses == 0.0)])})
            estimates[method][visits_string].update({'EMCI': np.mean(dpis[np.where(diagnoses == 0.25)])})
            estimates[method][visits_string].update({'LMCI': np.mean(dpis[np.where(diagnoses == 0.75)])})
            estimates[method][visits_string].update({'AD': np.mean(dpis[np.where(diagnoses == 1.0)])})

    for method in methods:
        print log.INFO, 'Results for {0}'.format(method)
        for diagnosis in ['CN', 'EMCI', 'LMCI', 'AD']:
            print log.RESULT, '{0: <4}:   {1:.2f} {2:.2f} | {3:.2f}  '.format(
                              diagnosis,
                              estimates[method]['m12'][diagnosis] - estimates[method]['bl'][diagnosis],
                              estimates[method]['m24'][diagnosis] - estimates[method]['m12'][diagnosis],
                              estimates[method]['m12_m24'][diagnosis] - estimates[method]['bl_m12'][diagnosis])

if __name__ == '__main__':
    main()
