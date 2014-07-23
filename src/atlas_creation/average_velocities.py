#! /usr/bin/env python
# print __doc__
import os.path
import argparse
from subprocess import call
from src.common import adni_tools as adni

EXEC_AVERAGE = 'ffdaverage'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    parser.add_argument('diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...')
    parser.add_argument('-t', '--trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic')
    a = parser.parse_args()

    dof_folder_adni1 = os.path.join(adni.data_folder, 'ADNI1', 'MNI152_' + a.trans + '_10mm_followup_to_baseline/dof')
    dof_folder_adni2 = os.path.join(adni.data_folder, 'ADNI2', 'MNI152_' + a.trans + '_10mm_followup_to_baseline/dof')

    velocities = adni.get_baseline_transformations(dof_folder_adni1, dof_folder_adni2, a.viscode, a.diagnosis)
    atlas_folder = os.path.join(adni.project_folder, 'atlas/model_0')
    out_average_velo = os.path.join(atlas_folder, 'velo_' + a.trans + '_' + a.viscode + '_' + a.diagnosis + '.dof.gz')

    print 'Found ' + str(len(velocities)) + ' velocities in total for viscode ' + a.viscode + '...'

    if os.path.isfile(out_average_velo):
        print 'Image already exists: ' + out_average_velo
    else:
        print '--------------------'
        print 'Output velocities: ' + out_average_velo

        call([EXEC_AVERAGE, out_average_velo] + velocities.tolist())


if __name__ == '__main__':
    main()
