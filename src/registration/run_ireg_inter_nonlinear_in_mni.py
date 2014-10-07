#! /usr/bin/env python2.7
import argparse
import os.path
import joblib as jl
from common import log as log
from common import adni_tools as adni
from common import atlas_tools as at
from registration import ireg_nonlinear


def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    # parser.add_argument('diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...')
    parser.add_argument('trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('state', type=float, help='the state for which relevant images should be registered')
    parser.add_argument('-n', '--nr_threads', type=int, default=1)
    # parser.add_argument('-i', '--iteration', type=int, default=1)
    parser.add_argument('-r', '--required_subjects', type=int, default=20)
    parser.add_argument('-s', '--spacing', type=str, default='10')
    parser.add_argument('-a', '--age_regression', action='store_true', default=False, help='use age regression')
    parser.add_argument('--save_image', action='store_true', help='save the warped image')
    a = parser.parse_args()

    global ireg_params
    ireg_params = os.path.join(adni.param_folder, 'params-ireg-' + a.trans + '-' + a.spacing + 'mm.txt')

    # atlas_folder = os.path.join(adni.project_folder, 'atlas')
    # datafile = os.path.join(atlas_folder, 'data_' + a.trans + '_' + a.viscode + '_' + a.diagnosis + '.csv')
    datafile = os.path.join(adni.project_folder, 'lists', 'dpis_for_atlas.csv')

    # rids, _, _, states, images = at.read_datafile(datafile, a.diagnosis, age_regression=a.age_regression)
    rids, states, images = at.read_datafile(datafile)

    _, _, indices = at.adaptive_kernel_regression(states, a.state, required_subjects=a.required_subjects,
                                                  sigma_min=1.0, sigma_max=50.0, sigma_delta=1.0)

    global selected_rids
    global selected_images
    selected_rids = rids[indices]
    selected_images = images[indices]

    global output_folder_img
    global output_folder_dof
    output_folder = adni.make_dir(adni.data_folder, 'ADNI/MNI152_intra_' + a.trans + '_' + a.spacing + 'mm')
    output_folder_dof = adni.make_dir(output_folder, 'dof')
    output_folder_img = adni.make_dir(output_folder, 'images')

    print log.RESULT, 'Found', len(selected_images), 'relevant images for state', a.state, '...'
    for j in range(len(selected_images)):
        jl.Parallel(n_jobs=a.nr_threads)(
            jl.delayed(run)(i, j, save_image=a.save_image) for i in range(len(selected_images)))


def run(index_targ, index_srce, save_image=False):
    target = selected_images[index_targ]
    source = selected_images[index_srce]
    target_rid = selected_rids[index_targ]
    source_rid = selected_rids[index_srce]
    out_basename = str(source_rid) + '_to_' + str(target_rid)

    out_dof = os.path.join(output_folder_dof, out_basename + '.dof.gz')
    if save_image:
        out_warped = os.path.join(output_folder_img, out_basename + '.nii.gz')
    else:
        out_warped = None

    ireg_nonlinear.run(source, target, None, out_dof, ireg_params, out_warped)


if __name__ == '__main__':
    main()
