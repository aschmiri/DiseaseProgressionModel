#! /usr/bin/env python2.7
import argparse
import os.path
from subprocess import call
from common import log as log
from common import adni_tools as adni

EXEC_NORM = 'asr_normalize_percentiles_apply'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('study')
    parser.add_argument('-a', '--algorithm', type=str, default='MNI152_linear')
    a = parser.parse_args()

    percentiles = '/vol/biomedic/users/reg09/posdoc/manifold_alignment/percentiles_10.csv'

    reg_folder = os.path.join(adni.data_folder, a.study, a.algorithm)
    image_folder = os.path.join(reg_folder, 'images')
    if not os.path.isdir(image_folder):
        print log.ERROR, 'Input folder does not exist.'
        exit()

    norm_folder = adni.make_dir(reg_folder, 'images_normalised')

    for infile in os.listdir(image_folder):
        if infile.endswith('.nii.gz'):
            image_in = os.path.join(image_folder, infile)
            image_out = os.path.join(norm_folder, infile)

            if os.path.exists(image_out):
                print log.SKIP, 'Image already exists: ', image_out
            else:
                print log.INFO, '--------------------'
                print log.INFO, 'Starting: normalize_percentiles_apply'
                print log.INFO, 'Input:  ' + image_in
                print log.INFO, 'Output: ' + image_out
                call([EXEC_NORM, percentiles, image_in, norm_folder + '/'])


if __name__ == '__main__':
    main()
