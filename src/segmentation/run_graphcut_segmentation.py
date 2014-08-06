#! /usr/bin/env python
import os.path
import argparse
import joblib as jl
from subprocess import call
from src.common import adni_tools as adni

EXEC_GRAPH = '/vol/biomedic/users/rw1008/linux_64/bin/atrophy_graphcut4D_v221'
EXEC_TRANS = 'transformation'
EXEC_THRESH = 'threshold'
EXEC_NORM = 'normalize'
EXEC_SCALE = 'imiImageShiftScale'


def main():
    parser = argparse.ArgumentParser(description='Perform the LEAP graphcut segmentation.')
    parser.add_argument('study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO')
    parser.add_argument('-n', '--nr_threads', dest='nr_threads', type=int, default=1)
    parser.add_argument('-l', '--lambda', dest='lmbda', type=float, default=15)
    parser.add_argument('-a', '--alpha', dest='alpha', type=float, default=0.2)
    parser.add_argument('-m', '--min_percentage', dest='min_percentage', type=float, default=0.84)
    parser.add_argument('-r', '--rid', dest='rid', type=int, default=None)
    args = parser.parse_args()

    pairs_collection = {}
    file_collection = adni.get_baseline_and_followups_as_collection()
    for rid in file_collection:
        if args.rid is None or rid == args.rid:
            bl_study = file_collection[rid]['bl']['study']
            if bl_study == args.study:
                bl_filename = file_collection[rid]['bl']['filename']
                bl_filepath = os.path.join(adni.data_folder, bl_study, 'native/images-1', bl_filename)
                fu_filepaths = []
                for viscode in file_collection[rid]:
                    if viscode != 'bl':
                        fu_study = file_collection[rid][viscode]['study']
                        fu_filename = file_collection[rid][viscode]['filename']
                        fu_filepath = os.path.join(adni.data_folder, fu_study, 'baseline_linear/images-1', fu_filename)
                        fu_filepaths.append(fu_filepath)
                pairs_collection.update({bl_filepath: fu_filepaths})

    print 'Found', len(pairs_collection), 'pairs...'
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(process_image_pair)(args, baseline, followups) for baseline, followups in pairs_collection.items())


def process_image_pair(args, baseline, followups):
    print '======================='
    print 'Starting: 4D graph cuts'
    print 'baseline: ', baseline
    print 'followup: ', followups

    # Define folders
    bl_study = adni.detect_study(baseline)
    out_dir = adni.make_dir(adni.data_folder, bl_study, 'native/seg_138regions_graphcut')
    image_dir = adni.make_dir(out_dir, 'images')

    # Get baseline name
    bl_name = os.path.basename(baseline).replace('.nii.gz', '')

    # Normalise followup images to baseline
    followups_normalised = []
    for followup in followups:
        fu_name = os.path.basename(followup).replace('.nii.gz', '')

        if not os.path.isfile(followup):
            print 'WARNING: warped follow up image not found:', followup
        else:
            fu_normalised = os.path.join(image_dir, fu_name + '_normalised.nii.gz')
            if not os.path.isfile(fu_normalised):
                print '---------------------'
                print 'Normalise follow-up image to baseline intensities'
                call([EXEC_NORM, baseline, followup, fu_normalised, '-Tp', '-1', '-Sp', '-1'])
            followups_normalised.append(fu_normalised)

    # Probabilistic atlases for tissue classes
    gm_mask = os.path.join(adni.mni_folder, 'mni_icbm152_gm_tal_nlin_asym_09a_scaled.nii')
    wm_mask = os.path.join(adni.mni_folder, 'mni_icbm152_wm_tal_nlin_asym_09a_scaled.nii')
    csf_mask = os.path.join(adni.mni_folder, 'mni_icbm152_csf_tal_nlin_asym_09a_scaled.nii')

    # Transform tissue priors to baseline space
    mni_dof_dir = os.path.join(adni.data_folder, bl_study, 'MNI152_linear', 'dof')
    mni_dof = os.path.join(mni_dof_dir, bl_name + '.dof.gz')
    if not os.path.isfile(mni_dof):
        print 'ERROR: transformation to MNI space not found:', mni_dof
        return

    print '---------------------'
    print 'Propagate tissue class priors for graph-cut optimisation'
    print 'MNI transformation: ', mni_dof
    gm_mask_transformed = os.path.join(image_dir, bl_name + '_gm.nii.gz')
    wm_mask_transformed = os.path.join(image_dir, bl_name + '_wm.nii.gz')
    csf_mask_transformed = os.path.join(image_dir, bl_name + '_csf.nii.gz')
    call([EXEC_TRANS, gm_mask, gm_mask_transformed, '-dofin', mni_dof, '-target', baseline, '-invert'])
    call([EXEC_TRANS, wm_mask, wm_mask_transformed, '-dofin', mni_dof, '-target', baseline, '-invert'])
    call([EXEC_TRANS, csf_mask, csf_mask_transformed, '-dofin', mni_dof, '-target', baseline, '-invert'])

    # Check atlas folder
    atlas_folder = os.path.join(adni.data_folder, bl_study, 'native', 'seg_138regions_probs', bl_name)
    if not os.path.isdir(atlas_folder):
        print 'ERROR: Atlas folder not found:', atlas_folder
        return

    # Process all structures individually
    for structure in adni.volume_names:
        # Get atlas for structure
        index = adni.volume_names.index(structure) + 1
        atlas = os.path.join(image_dir, bl_name + '_posteriors_{0}.nii.gz'.format(index))
        if not os.path.isfile(atlas):
            print '---------------------'
            print 'Scaling atlas to [0..100]'
            atlas_unscaled = os.path.join(atlas_folder, 'posteriors_{0}.nii.gz'.format(index))
            call([EXEC_SCALE, '-I', atlas_unscaled, '-O', atlas, '-a', '0', '-m', '100'])

        # Perform graph cut
        print '---------------------'
        print 'Perform graph-cut optimisation to obtain final segmentation'
        call([EXEC_GRAPH, str(len(followups_normalised) + 1), baseline] + followups_normalised +
             ['1', atlas, wm_mask_transformed, gm_mask_transformed, csf_mask_transformed,
              str(args.lmbda), str(args.alpha), out_dir, '-minPercentage', str(args.min_percentage)])

        # Rename output segmentations
        postfix = '_' + structure.replace(' ', '_')
        bl_out = os.path.join(out_dir, bl_name + postfix + '.nii.gz')
        call(['mv', os.path.join(out_dir, '{0}_0.nii.gz'.format(bl_name)), bl_out])
        for followup in followups_normalised:
            fu_name = os.path.basename(followup).replace('_normalised.nii.gz', '')
            fu_out = os.path.join(out_dir, fu_name + postfix + '.nii.gz')
            fu_index = followups_normalised.index(followup) + 1
            call(['mv', os.path.join(out_dir, '{0}_{1}.nii.gz'.format(bl_name, fu_index)), fu_out])


if __name__ == '__main__':
    main()
