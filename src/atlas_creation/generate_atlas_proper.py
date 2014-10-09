#! /usr/bin/env python2.7
import os.path
import argparse
import csv
from subprocess import call
from common import log as log
from common import adni_tools as adni
from common import atlas_tools as at

EXEC_ATLAS = 'atlas'
EXEC_FFDAVERAGE = 'ffdaverage'
EXEC_TRANSFORM = 'transformation'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('state', type=float, help='the state for which relevant images should be registered')
    parser.add_argument('-r', '--required_subjects', type=int, default=20)
    parser.add_argument('-t', '--trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('-s', '--spacing', type=str, default='10')
    parser.add_argument('-a', '--age_regression', action='store_true', default=False, help='use age regression')
    a = parser.parse_args()

    image_method_folder = 'MNI152_linear'
    image_type_folder = 'images_normalised'

    dof_folder = os.path.join(adni.data_folder, 'ADNI', 'MNI152_intra_' + a.trans + '_' + a.spacing + 'mm', 'dof')

    atlas_folder = os.path.join(adni.project_folder, 'atlas/all_no_cn')
    atlas_folder_temp = adni.make_dir(atlas_folder, 'temp')
    datafile = os.path.join(adni.project_folder, 'lists', 'dpis_for_atlas.csv')

    # Find sigma and corresponding images
    rids, states, images = at.read_datafile(datafile)
    _, weights, indices = at.adaptive_kernel_regression(states, a.state, required_subjects=a.required_subjects,
                                                        sigma_min=1.0, sigma_max=50.0, sigma_delta=1.0)

    selected_rids = rids[indices]
    selected_images = images[indices]
    selected_weights = weights[indices]

    # Print data file for IRTK image averaging
    atlas_base = 'atlas_' + a.trans + '_' + str(a.state)

    if a.age_regression:
        atlas_base = atlas_base + '_agereg'
    data_file_images = os.path.join(atlas_folder, atlas_base + '_images.txt')

    with open(data_file_images, 'wb') as csvfile:
        image_writer = csv.writer(csvfile, delimiter=' ')
        for i in range(len(selected_images)):
            target = at.find_file(selected_images[i], image_method_folder, image_type_folder)
            target_rid = selected_rids[i]

            # Define the base name of output files
            temp_base = str(target_rid) + '_' + str(a.state)
            if a.age_regression:
                temp_base = temp_base + '_agereg'

            # Print data file for IRTK ffd averaging
            data_file_dofs = os.path.join(atlas_folder_temp, temp_base + '_dofs.txt')
            with open(data_file_dofs, 'wb') as csvfile:
                dof_writer = csv.writer(csvfile, delimiter=' ')
                for j in range(len(selected_images)):
                    if i != j:
                        source_rid = selected_rids[j]
                        dof_basename = str(source_rid) + '_to_' + str(target_rid) + '.dof.gz'
                        dof = os.path.join(dof_folder, dof_basename)
                        dof_writer.writerow([dof, selected_weights[j]])

            # Call IRTK 'ffdaverage'
            out_average_dof = os.path.join(atlas_folder_temp, temp_base + '_average_ffd.dof.gz')
            if os.path.exists(out_average_dof):
                print log.SKIP, 'Transformation {0} already exists.'.format(out_average_dof)
            else:
                print log.INFO, '--------------------'
                print log.INFO, 'Starting: ffdaverage'
                print log.INFO, 'Data:   ' + data_file_dofs
                print log.INFO, 'Output: ' + out_average_dof
                call([EXEC_FFDAVERAGE, out_average_dof, '-dofnames', data_file_dofs,
                      '-identity', str(selected_weights[i])])

            # Transform average image
            out_image_transformed = os.path.join(atlas_folder_temp, temp_base + '_image_transformed.nii.gz')
            if os.path.exists(out_image_transformed):
                print log.SKIP, 'Transformed image {0} already exists.'.format(out_average_dof)
            else:
                print log.INFO, '--------------------'
                print log.INFO, 'Starting: transformation'
                print log.INFO, 'Input:  ' + target
                print log.INFO, 'DOF:    ' + out_average_dof
                print log.INFO, 'Output: ' + out_image_transformed
                call([EXEC_TRANSFORM, target, out_image_transformed,
                      '-dofin', out_average_dof, '-invert'])

            # Write transformed image to image file
            image_writer.writerow([out_image_transformed, selected_weights[i]])

    # Call IRTK 'atlas'
    out_average_image = os.path.join(atlas_folder, atlas_base + '_average_image.nii.gz')
    if os.path.exists(out_average_image):
        print log.SKIP, 'Atlas {0} already exists.'.format(out_average_dof)
    else:
        print log.INFO, '--------------------'
        print log.INFO, 'Starting: ffdaverage'
        print log.INFO, 'Data:   ' + data_file_images
        print log.INFO, 'Output: ' + out_average_image
        call([EXEC_ATLAS, out_average_image, '-imagenames', data_file_images])

    print log.RESULT, 'Atlas generated successfully.'


if __name__ == '__main__':
    main()
