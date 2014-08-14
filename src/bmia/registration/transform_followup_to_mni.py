#! /usr/bin/env python2.7
import argparse
import os.path
import joblib as jl
from subprocess import call
from bmia.common import log as log
from bmia.common import adni_tools as adni


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO')
    parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    parser.add_argument('trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('-n', '--nr_threads', dest='nr_threads', type=int, default=1)
    parser.add_argument('-s', '--spacing', dest='sx', type=str, default='10')
    a = parser.parse_args()

    global mni_linear_folder_dof
    global mni_nonlin_folder_dof
    data_folder = os.path.join(adni.data_folder, a.study)

    mni_linear_folder = os.path.join(data_folder, 'MNI152_linear')
    mni_linear_folder_dof = os.path.join(mni_linear_folder, 'dof')

    mni_nonlin_folder = os.path.join(data_folder, 'MNI152_' + a.trans + '_' + a.sx + 'mm_after_linear')
    mni_nonlin_folder_dof = os.path.join(mni_nonlin_folder, 'dof')

    global out_folder_linear_img
    global out_folder_nonlin_img
    out_folder_linear = adni.make_dir(data_folder, 'MNI152_linear_via_baseline')
    out_folder_linear_img = adni.make_dir(out_folder_linear, 'images')
    out_folder_nonlin_img = adni.make_dir(mni_nonlin_folder, 'images_' + a.viscode)

    global baseline_files
    global followup_files
    baseline_folder = os.path.join(data_folder, 'native/images')
    followup_folder = os.path.join(data_folder, 'baseline_linear/images')
    baseline_files, followup_files = adni.get_baseline_and_followup(baseline_folder, followup_folder, a.study, a.viscode)

    print log.RESULT, 'Found', len(followup_files), 'images...'
    jl.Parallel(n_jobs=a.nr_threads)(jl.delayed(run, a.study)(i) for i in range(len(followup_files)))


def run(index, study):
    baseline = baseline_files[index]
    followup = followup_files[index]
    baseline_base = os.path.basename(baseline)
    followup_base = os.path.basename(followup)
    aff = os.path.join(mni_linear_folder_dof, baseline_base.replace('.nii.gz', '.dof.gz'))
    dof = os.path.join(mni_nonlin_folder_dof, baseline_base.replace('.nii.gz', '.dof.gz'))

    out_image_affine = os.path.join(out_folder_linear_img, followup_base)
    out_image_nonlin = os.path.join(out_folder_nonlin_img, followup_base)

    if os.path.isfile(out_image_affine):
        print log.SKIP, 'Image already exists: ' + out_image_nonlin
    else:
        print log.INFO, '--------------------'
        print log.INFO, 'In image: ', followup
        print log.INFO, 'Affine:   ', aff
        print log.INFO, 'Out image:', out_image_affine
        call(['transformation', followup, out_image_affine, '-dofin', aff, '-target', adni.mni_atlas, '-cspline', '-matchInputType', '-Sp', '0'])

    if os.path.isfile(out_image_nonlin):
        print log.SKIP, 'Image already exists: ' + out_image_nonlin
    else:
        print log.INFO, '--------------------'
        print log.INFO, 'In image: ', out_image_affine
        print log.INFO, 'Nonlinear:', dof
        print log.INFO, 'Out image:', out_image_nonlin
        call(['transformation', out_image_affine, out_image_nonlin, '-dofin', dof, '-target', adni.mni_atlas, '-cspline', '-matchInputType', '-Sp', '0'])


if __name__ == '__main__':
    main()
