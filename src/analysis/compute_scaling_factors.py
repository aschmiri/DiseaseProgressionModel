#! /usr/bin/env python
# print __doc__
import argparse
import csv
import os.remove
from subprocess import check_output
from src.common import adni_tools as adni

EXEC_FACTORS = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/affineScaling'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, default='/vol/medic01/users/aschmidt/projects/Data/ADNI/data/ADNIX/MNI152_linear/dof', help='the path with the dof files with the affine transformations to MNI space')
    parser.add_argument('--output', type=str, default='./scaling_factors.csv', help='the output file with the scaling factors')
    a = parser.parse_args()

    os.remove(a.output)
    compute_for_study('ADNI1', a.folder, a.output)
    compute_for_study('ADNIGO', a.folder, a.output)
    compute_for_study('ADNI2', a.folder, a.output)


def compute_for_study(study, folder, output):
    dof_folder = folder.replace('ADNIX', study)
    files, rids, _, _, _ = adni.read_list_all_data(dof_folder, diagnosis='ALL', study=study, viscode='bl')

    print 'Found', str(len(files)), 'images for', study, '...'
    with open(output, 'a') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['RID', 'Factor', 'Filename'])
        for rid, dof in zip(rids, files):
            factor = float(check_output([EXEC_FACTORS, dof]))
            writer.writerow([rid, factor, dof])


if __name__ == '__main__':
    main()
