#! /usr/bin/env python2.7
import argparse
import os.path
from subprocess import call
from bmia.common import log as log
from bmia.common import adni_tools as adni

EXEC_RVIEW = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/rview'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO')
    parser.add_argument('trans', type=str, help='the transformation model, e.g. linear, ffd, svffd, sym, or ic')
    parser.add_argument('-d', '--dof', action='store_true', default=False, help='show the dof')
    parser.add_argument('-r', '--rid', type=int, default=None)
    parser.add_argument('-s', '--spacing', dest='sx', type=str, default='10')
    a = parser.parse_args()

    data_folder = os.path.join(adni.data_folder, a.study)
    if a.trans == 'linear':
        if a.dof:
            baseline_folder = os.path.join(data_folder, 'native/images')
        else:
            baseline_folder = os.path.join(data_folder, 'MNI152_linear/images')
        dof_folder = os.path.join(data_folder, 'MNI152_linear/dof')
    else:
        if a.dof:
            baseline_folder = os.path.join(data_folder, 'MNI152_linear/images')
        else:
            baseline_folder = os.path.join(data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear/images')
        dof_folder = os.path.join(data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear/dof')

    print dof_folder
    print baseline_folder
    baseline_files = adni.get_baseline(baseline_folder, a.study)
    baseline_files, dof_files = adni.find_images_with_dof(baseline_files, dof_folder)

    print log.RESULT, 'Found ' + str(len(baseline_files)) + ' images:'
    for i in range(len(baseline_files)):
        source = baseline_files[i]
        dof = dof_files[i]
        if a.rid is None or adni.detect_rid(source) == a.rid:
            print log.INFO, '--------------------'
            print log.INFO, 'Source: ' + source
            print log.INFO, 'DOF:    ' + dof

            if a.dof:
                call([EXEC_RVIEW, adni.mni_atlas, source, dof, '-res', '1.5', '-mix'])
            else:
                call([EXEC_RVIEW, adni.mni_atlas, source, '-res', '1.5', '-mix'])


if __name__ == '__main__':
    main()
