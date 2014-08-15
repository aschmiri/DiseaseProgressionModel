#! /usr/bin/env python2.7
import os.path
import argparse
import joblib as jl
from subprocess import call
import numpy as np
from common import log as log
from common import adni_tools as adni

EXEC_GRAPH = '/vol/biomedic/users/rw1008/linux_64/bin/atrophy_graphcut4D_v221'
EXEC_TRANS = 'transformation'
EXEC_THRESH = 'threshold'
EXEC_NORM = 'normalize'
EXEC_SCALE = 'imiImageShiftScale'


def main():
    parser = argparse.ArgumentParser(description='Perform the LEAP graphcut segmentation.')
    parser.add_argument('-n', '--nr_threads', dest='nr_threads', type=int, default=1)
    parser.add_argument('-l', '--lambda', dest='lmbda', type=float, default=15)
    parser.add_argument('-a', '--alpha', dest='alpha', type=float, default=0.2)
    parser.add_argument('-m', '--min_percentage', dest='min_percentage', type=float, default=0.84)
    parser.add_argument('--study', type=str, default='ALL', help='the study, should be ADNI1, ADNI2, or ADNIGO')
    parser.add_argument('--min_scans', type=int, default=1, help='the minimal number of scans per subject')
    parser.add_argument('--rid', type=int, default=None)
    parser.add_argument('--rid_mod', type=int, default=None)
    parser.add_argument('--all_structures', action='store_true', default=False, help='segment all 138 structures')
    args = parser.parse_args()

    file_collection = adni.get_baseline_and_followups_as_collection()

    print log.INFO, 'Assembling list of subjects...'
    pairs_collection = {}
    for rid in file_collection:
        if (args.rid is None and args.rid_mod is None) or (args.rid_mod is not None and np.mod(rid, 5) == args.rid_mod) or rid == args.rid:
            bl_study = file_collection[rid]['bl']['study']
            if args.study == 'ALL' or bl_study == args.study:
                # Collect baseline data
                bl_filename = file_collection[rid]['bl']['filename']
                bl_filepath = os.path.join(adni.data_folder, bl_study, 'native/images-1', bl_filename)

                # Collect followup data
                fu_filepaths = []
                fu_scantimes = []
                for viscode in file_collection[rid]:
                    if viscode != 'bl':
                        fu_study = file_collection[rid][viscode]['study']
                        fu_filename = file_collection[rid][viscode]['filename']
                        fu_filepath = os.path.join(adni.data_folder, fu_study, 'baseline_linear/images-1', fu_filename)
                        fu_filepaths.append(fu_filepath)
                        fu_scantimes.append(int(viscode[1:]))

                if len(fu_filepaths) >= args.min_scans:
                    # Sort followup data
                    fu_filepaths = np.array(fu_filepaths)
                    fu_scantimes = np.array(fu_scantimes)

                    indices = np.argsort(fu_scantimes)
                    fu_filepaths = fu_filepaths[indices]
                    fu_scantimes = fu_scantimes[indices]

                    # Add to pairs collection
                    pairs_collection.update({bl_filepath: fu_filepaths})

    print log.RESULT, 'Found', len(pairs_collection), 'pairs...'
    jl.Parallel(n_jobs=args.nr_threads)(jl.delayed(process_image_pairs)(args, baseline, followups) for baseline, followups in pairs_collection.items())


def process_image_pairs(args, baseline, followups):
    print log.INFO, '======================='
    print log.INFO, 'STARTING 4D graph cuts'
    print log.INFO, 'baseline: ', os.path.basename(baseline)
    for i, followup in enumerate(followups):
        print log.INFO, 'followup {0}: {1}'.format(i, os.path.basename(followup))

    # Determine structures to segment
    if args.all_structures:
        structures = adni.volume_names
    else:
        structures = adni.volume_names_essential

    # Define folders
    bl_study = adni.detect_study(baseline)
    out_dir = adni.make_dir(adni.data_folder, bl_study, 'native/seg_138regions_graphcut')
    temp_dir = adni.make_dir(out_dir, 'temp')

    # Get baseline name
    bl_name = os.path.basename(baseline).replace('.nii.gz', '')

    # Normalise followup images to baseline
    followups_normalised = []
    for followup in followups:
        if not os.path.isfile(followup):
            print log.WARNING, 'No warped follow up image found:', followup
        else:
            fu_name = os.path.basename(followup).replace('.nii.gz', '')
            fu_normalised = os.path.join(temp_dir, fu_name + '_normalised.nii.gz')
            if os.path.isfile(fu_normalised):
                print log.SKIP, 'Warped followup available for {0}'.format(bl_name)
            else:
                if not os.path.isfile(fu_normalised):
                    print log.INFO, 'Normalising follow-up image to baseline intensities...'
                    call([EXEC_NORM, baseline, followup, fu_normalised, '-Tp', '-1', '-Sp', '-1'])
            followups_normalised.append(fu_normalised)

    # Probabilistic atlases for tissue classes
    gm_mask_transformed = os.path.join(temp_dir, bl_name + '_gm.nii.gz')
    wm_mask_transformed = os.path.join(temp_dir, bl_name + '_wm.nii.gz')
    csf_mask_transformed = os.path.join(temp_dir, bl_name + '_csf.nii.gz')
    if files_present([gm_mask_transformed, wm_mask_transformed, csf_mask_transformed]):
        print log.SKIP, 'Tissue class priors available for {0}'.format(bl_name)
    else:
        print log.INFO, 'Propagating tissue class priors for graph-cut optimisation...'
        gm_mask = os.path.join(adni.mni_folder, 'mni_icbm152_gm_tal_nlin_asym_09a_scaled.nii')
        wm_mask = os.path.join(adni.mni_folder, 'mni_icbm152_wm_tal_nlin_asym_09a_scaled.nii')
        csf_mask = os.path.join(adni.mni_folder, 'mni_icbm152_csf_tal_nlin_asym_09a_scaled.nii')

        # Transform tissue priors to baseline space
        mni_dof_dir = os.path.join(adni.data_folder, bl_study, 'MNI152_linear', 'dof')
        mni_dof = os.path.join(mni_dof_dir, bl_name + '.dof.gz')
        if not os.path.isfile(mni_dof):
            print log.ERROR, 'Transformation to MNI space not found:', mni_dof
            return

        call([EXEC_TRANS, gm_mask, gm_mask_transformed, '-dofin', mni_dof, '-target', baseline, '-invert'])
        call([EXEC_TRANS, wm_mask, wm_mask_transformed, '-dofin', mni_dof, '-target', baseline, '-invert'])
        call([EXEC_TRANS, csf_mask, csf_mask_transformed, '-dofin', mni_dof, '-target', baseline, '-invert'])

    # Check atlas folder
    atlas_folder = os.path.join(adni.data_folder, bl_study, 'native', 'seg_138regions_probs', bl_name)
    if not os.path.isdir(atlas_folder):
        print log.ERROR, 'Atlas folder not found:', atlas_folder
        return

    # Process all structures individually
    for structure in structures:
        output_names = get_output_names(out_dir, bl_name, followups_normalised, structure)

        if files_present(output_names):
            print log.SKIP, 'Segmentations of {0} present for {1}'.format(structure, bl_name)
        else:
            print log.INFO, '---------------------'
            print log.INFO, 'Segmenting {0}...'.format(structure)

            # Get atlas for structure
            posterior_index = adni.volume_names.index(structure) + 1
            atlas = os.path.join(temp_dir, bl_name + '_posteriors_{0}.nii.gz'.format(posterior_index))
            if not os.path.isfile(atlas):
                print log.INFO, 'Scaling atlas to [0..100]...'
                atlas_unscaled = os.path.join(atlas_folder, 'posteriors_{0}.nii.gz'.format(posterior_index))
                call([EXEC_SCALE, '-I', atlas_unscaled, '-O', atlas, '-a', '0', '-m', '100'])

            # Perform graph cut
            print log.INFO, 'Performing graph-cut optimisation...'
            call([EXEC_GRAPH, str(len(followups_normalised) + 1), baseline] + followups_normalised +
                 ['1', atlas, wm_mask_transformed, gm_mask_transformed, csf_mask_transformed,
                  str(args.lmbda), str(args.alpha), out_dir, '-minPercentage', str(args.min_percentage)])

            # Rename output segmentations
            print log.INFO, 'Renaming files...'
            for i, output_name in enumerate(output_names):
                call(['mv', os.path.join(out_dir, '{0}_{1}.nii.gz'.format(bl_name, i)), output_name])

    # Clean temporal files
    # call('rm ' + os.path.join(temp_dir, bl_name[:15] + '*'), shell=True)

    print log.INFO, 'Finished processing of {0}'.format(bl_name)
    print log.INFO, '======================='


def get_output_names(out_dir, bl_name, fu_names, structure):
    postfix = '_' + structure.replace(' ', '_')
    output_names = [os.path.join(out_dir, bl_name + postfix + '.nii.gz')]
    for followup in fu_names:
        fu_name = os.path.basename(followup).replace('_normalised', '').replace('.nii.gz', '')
        output_names.append(os.path.join(out_dir, fu_name + postfix + '.nii.gz'))
    return output_names


def files_present(filenames):
    all_found = True
    for filename in filenames:
        all_found = all_found and os.path.isfile(filename)
    return all_found

if __name__ == '__main__':
    main()
