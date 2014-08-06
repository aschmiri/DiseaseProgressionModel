#! /usr/bin/env python
# print __doc__
import argparse
import os.path
import joblib as jl
from subprocess import call
from src.common import adni_tools as adni

EXEC_MASK = 'padding'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO')
    parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    parser.add_argument('-n', '--nr_threads', dest='nr_threads', type=int, default=1)
    args = parser.parse_args()

    global mask_folder
    data_folder = os.path.join(adni.data_folder, args.study)
    mask_folder = os.path.join(data_folder, 'native/masks_brain')

    global baseline_files
    global followup_files
    baseline_folder = os.path.join(data_folder, 'native/images_unstripped')
    followup_folder = os.path.join(data_folder, 'baseline_linear/images_unstripped')
    baseline_files, followup_files = adni.get_baseline_and_followup(baseline_folder, followup_folder, args.study, args.viscode)

    global output_folder
    output_folder = adni.make_dir(data_folder, 'baseline_linear/images')

    print 'Found ' + str(len(baseline_files)) + ' images...'
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(run)(i) for i in range(len(baseline_files)))


def run(index):
    for i in range(len(baseline_files)):
        baseline = baseline_files[i]
        baseline_base = os.path.basename(baseline)
        followup = followup_files[i]
        followup_base = os.path.basename(followup)

        mask = os.path.join(mask_folder, baseline_base)
        mask = adni.find_file(mask)
        if mask is not None:
            out = os.path.join(output_folder, followup_base)

            if os.path.isfile(out):
                print 'Image already exists: ' + out
            elif mask is None or not os.path.isfile(mask):
                print 'No mask found for: ' + out
            else:
                print '--------------------'
                print 'Image:  ' + followup
                print 'Mask:   ' + mask
                print 'Output: ' + out
                call([EXEC_MASK, followup, mask, out, '1', '0', '-invert'])


if __name__ == '__main__':
    main()
