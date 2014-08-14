#! /usr/bin/env python2.7
import os.path
import argparse
import csv
from subprocess import check_output
from bmia.common import adni_tools as adni

EXEC_ESTIMATE = 'stateestimation'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    parser.add_argument('diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...')
    parser.add_argument('-i', '--iteration', dest='iteration', type=int, default=0)
    parser.add_argument('-t', '--trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('--state_min', type=float, default=0)
    parser.add_argument('--state_max', type=float, default=15)
    parser.add_argument('--state_stepsize', type=float, default=0.25)
    parser.add_argument('--postfix', type=str, default='')
    a = parser.parse_args()

    mask_brain = os.path.join(adni.mni_folder, 'MNI152_T1_1mm_brainmask.nii')

    atlas_folder = os.path.join(adni.project_folder, 'atlas/model_' + str(a.iteration) + a.postfix)

    image_folder_adni1 = os.path.join(adni.data_folder, 'ADNI1', 'MNI152_linear/images')
    image_folder_adni2 = os.path.join(adni.data_folder, 'ADNI2', 'MNI152_linear/images')
    image_folder_adniG = os.path.join(adni.data_folder, 'ADNIGO', 'MNI152_linear/images')
    datafile = os.path.join(atlas_folder, 'data_' + a.trans + '_' + a.viscode + '_' + a.diagnosis + '.csv')
    model_prefix = os.path.join(atlas_folder, 'model_' + a.trans + '_' + a.viscode + '_' + a.diagnosis + '_p')

    files, rids, diags, ages, mmses = adni.get_all_data(
        image_folder_adni1, image_folder_adni2, image_folder_adniG, 'bl')

    if os.path.exists(datafile):
        print ' File ' + datafile + ' already exists.'
    else:
        with open(datafile, 'wb') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',')
            csv_writer.writerow(['RID', 'DX.bl', 'AGE', 'MMSE', 'DPI', 'FILE'])
            for i in range(len(files)):
                print ' Analysing ' + files[i]

                # Estimate disease state for image
                if os.path.exists(files[i]) and a.diagnosis in diags[i] or diags[i] in a.diagnosis:
                    virtual_nmi_brain = check_output([EXEC_ESTIMATE, files[i], '-mask', mask_brain,
                                                      '-min', str(a.state_min), '-max', str(a.state_max), '-stepsize', str(a.state_stepsize),
                                                      '-prefix', model_prefix, '-silent'])
                    csv_writer.writerow([rids[i], diags[i], ages[i], mmses[i], virtual_nmi_brain, files[i]])
                    print virtual_nmi_brain


if __name__ == '__main__':
    main()
