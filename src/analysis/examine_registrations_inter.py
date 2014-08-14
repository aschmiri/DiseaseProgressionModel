#! /usr/bin/env python
# print __doc__
import argparse
import os.path
from subprocess import call
from src.common import adni_tools as adni
from src.common import atlas_tools as at

EXEC_RVIEW = 'rview'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('state', type=float, help='the state for which relevant images should be registered')
    parser.add_argument('-t', '--target_rid', type=int, default=None)
    parser.add_argument('-s', '--spacing', dest='sx', type=str, default='10')
    parser.add_argument('-r', '--required_subjects', dest='required_subjects', type=int, default=20)
    a = parser.parse_args()

    datafile = os.path.join(adni.project_folder, 'atlas/model_1/data_m24_AD.csv')
    rids, _, _, states, images = at.read_datafile(datafile, 'AD')

    _, _, indices = at.adaptive_kernel_regression(states, a.state, required_subjects=a.required_subjects)

    selected_rids = rids[indices]
    selected_images = images[indices]

    dof_folder = os.path.join(adni.data_folder, 'ADNI/MNI152_intra_' + a.trans + '_' + a.sx + 'mm', 'dof')

    print adni.RESULT, 'Found ' + str(len(selected_images)) + ' relevant images for state ' + str(a.state) + '...'
    for i in range(len(selected_images)):
        for j in range(len(selected_images)):
            target = selected_images[i]
            source = selected_images[j]
            target_rid = selected_rids[i]
            source_rid = selected_rids[j]
            dof_basename = str(source_rid) + '_to_' + str(target_rid)

            dof = os.path.join(dof_folder, dof_basename + '.dof.gz')

            if os.path.exists(dof):
                if a.target_rid is None or a.target_rid == target_rid:
                    print adni.INFO, '--------------------'
                    print adni.INFO, 'Target: ' + target
                    print adni.INFO, 'Source: ' + source
                    print adni.INFO, 'DOF:    ' + dof

                    call([EXEC_RVIEW, target, source, dof, '-res', '1.5', '-mix'])


if __name__ == '__main__':
    main()
