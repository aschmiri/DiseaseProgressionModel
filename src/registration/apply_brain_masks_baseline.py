#! /usr/bin/env python2.7
import argparse
import os.path
import joblib as jl
from subprocess import call
from common import log as log
from common import adni_tools as adni


EXEC_MASK = 'padding'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO')
    parser.add_argument('-n', '--nr_threads', dest='nr_threads', type=int, default=1)
    parser.add_argument('-p', '--padding', dest='padding', type=int, default=0)
    args = parser.parse_args()

    global mask_folder
    data_folder = os.path.join(adni.data_folder, args.study)
    mask_folder = os.path.join(data_folder, 'native/masks_brain')

    global baseline_files
    baseline_folder = os.path.join(data_folder, 'native/images_unstripped')
    baseline_files = adni.get_baseline(baseline_folder, args.study)

    global output_folder
    output_folder = adni.make_dir(data_folder, 'native/images')

    print log.RESULT, 'Found', len(baseline_files), 'images...'
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(run)(args, i) for i in range(len(baseline_files)))


def run(args, index):
    image = baseline_files[index]
    image_base = os.path.basename(image)

    mask = os.path.join(mask_folder, image_base)
    mask = adni.find_alternative_file(mask)
    if mask is not None:
        out = os.path.join(output_folder, image_base)

        if os.path.isfile(out):
            print log.SKIP, 'Image already exists:', out
        elif not os.path.isfile(mask):
            print log.INFO, 'No mask found for: ', out
        else:
            print log.INFO, '--------------------'
            print log.INFO, 'Image: ', image
            print log.INFO, 'Mask:  ', mask
            print log.INFO, 'Output:', out

            call([EXEC_MASK, image, mask, out, '1', str(args.padding), '-invert'])


if __name__ == '__main__':
    main()
