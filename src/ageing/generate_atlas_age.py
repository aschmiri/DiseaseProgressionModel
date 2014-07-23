#! /usr/bin/env python
# print __doc__
import argparse
import os.path
import csv
from subprocess import call
from src.common import adni_tools as adni
from src.common import atlas_tools as at

EXEC_ATLAS = 'atlas'
EXEC_FFDAVERAGE = 'ffdaverage'
EXEC_TRANSFORM = 'transformation'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    parser.add_argument('diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...')
    parser.add_argument('age', type=float, help='the age for which relevant images should be registered')
    parser.add_argument('-r', '--required_subjects', type=int, default=20)
    parser.add_argument('-t', '--trans', type=str, default='svffd', help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('-s', '--spacing', type=str, default='10')
    a = parser.parse_args()

    image_method_folder = 'MNI152_linear'
    image_type_folder = 'images_normalised'

    data_folder = '/vol/biomedic/users/aschmidt/ADNI/data/ADNI'
    dof_folder = os.path.join(data_folder, 'MNI152_intra_' + a.trans + '_' + a.spacing + 'mm', 'dof')

    atlas_folder = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/ageing'
    atlas_folder_temp = adni.make_dir(atlas_folder, 'temp')
    # datafile = os.path.join( atlas_folder, 'data_' + a.viscode + '_' + a.diagnosis + '.csv' )

    data_folder = '/vol/biomedic/users/aschmidt/ADNI'
    image_folder_adni1 = os.path.join(data_folder, 'data/ADNI1/MNI152_linear/images_normalised')
    image_folder_adni2 = os.path.join(data_folder, 'data/ADNI2/MNI152_linear/images_normalised')
    image_folder_adniG = os.path.join(data_folder, 'data/ADNIGO/MNI152_linear/images_normalised')

    images, rids, _, ages, _ = adni.get_all_data(
            image_folder_adni1, image_folder_adni2, image_folder_adniG, a.viscode, diagnosis=a.diagnosis)

    # Find sigma and corresponding images
    _, weights, indices = at.adaptive_kernel_regression(ages, a.age, required_subjects=a.required_subjects)

    selected_rids = rids[indices]
    selected_images = images[indices]
    selected_weights = weights[indices]

    # Print data file for IRTK image averaging
    atlas_base = 'atlas_' + a.viscode + '_' + a.diagnosis + '_' + str(a.age)
    data_file_images = os.path.join(atlas_folder, atlas_base + '_images.txt')

    with open(data_file_images, 'wb') as csvfile:
        image_writer = csv.writer(csvfile, delimiter=' ')
        for i in range(len(selected_images)):
            target = at.find_file(selected_images[i], image_method_folder, image_type_folder)
            target_rid = selected_rids[i]

            # Define the base name of output files
            temp_base = str(target_rid) + '_' + str(a.age) + '_' + a.diagnosis + '_' + str(a.required_subjects)

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
            if not os.path.exists(out_average_dof):
                print '--------------------'
                print 'Starting: ffdaverage'
                print 'Data:   ' + data_file_dofs
                print 'Output: ' + out_average_dof
                call([EXEC_FFDAVERAGE, out_average_dof, '-dofnames', data_file_dofs,
                       '-identity', str(selected_weights[i])])

            # Transform average image
            out_image_transformed = os.path.join(atlas_folder_temp, temp_base + '_image_transformed.nii.gz')
            if not os.path.exists(out_image_transformed):
                print '--------------------'
                print 'Starting: transformation'
                print 'Input:  ' + target
                print 'DOF:    ' + out_average_dof
                print 'Output: ' + out_image_transformed
                call([EXEC_TRANSFORM, target, out_image_transformed,
                       '-dofin', out_average_dof, '-invert'])  # , '-cspline'

            # Write transformed image to image file
            image_writer.writerow([out_image_transformed, selected_weights[i]])

    # Call IRTK 'atlas'
    out_average_image = os.path.join(atlas_folder, atlas_base + '_average_image.nii.gz')
    if not os.path.exists(out_average_image):
        print '--------------------'
        print 'Starting: ffdaverage'
        print 'Data:   ' + data_file_images
        print 'Output: ' + out_average_image
        call([EXEC_ATLAS, out_average_image, '-imagenames', data_file_images])


if __name__ == '__main__':
    main()
