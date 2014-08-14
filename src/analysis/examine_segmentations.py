#! /usr/bin/env python
# print __doc__
import argparse
import os.path
import sqlite3
from subprocess import call
from src.common import adni_tools as adni

EXEC_RVIEW = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/rview'
EXEC_EXTRACT = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_extract_frame'
EXEC_REMOVE = 'rm'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('method', choices=['reg', 'long', 'cons', 'graph'])
    parser.add_argument('-r', '--rid', type=int, default=None)
    parser.add_argument('-s', '--structure_index', type=int, default=None)
    parser.add_argument('-v', '--viscode', type=str, default=None)
    parser.add_argument('--trans', dest='trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic (regbased only)')
    parser.add_argument('--spacing', dest='sx', type=str, default='5', help='the transformation spacing (regbased only)')
    parser.add_argument('--all_structures', action='store_true', default=False, help='segment all 138 structures (graphcuts only)')
    args = parser.parse_args()

    lut = '/vol/medic02/users/cl6311/Neuro_Atlas/config/lut.csv'

    # Determine structures to segment
    if args.method == 'graph' and not args.all_structures:
        structures = adni.volume_names_essential
    else:
        structures = adni.volume_names

    # Setup DB
    con = sqlite3.connect(os.path.join(adni.project_folder, 'lists', 'adni.db'))
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    # Get all RIDs
    cur.execute("SELECT rid FROM adnimerge WHERE viscode = 'bl'")
    rows = cur.fetchall()
    rids = [row['rid'] for row in rows]

    for rid in rids:
        if args.rid is not None and rid != args.rid:
            continue

        # Get all visits for RID
        cur.execute("SELECT iid, viscode, study, study_bl, diagnosis, filename, age, scandate \
                      FROM adnimerge WHERE rid = " + str(rid))
        visits = cur.fetchall()

        # Get sorted viscodes (only needed for consistent 4D segmentation)
        sorted_viscodes = get_sorted_viscodes(cur, rid)

        for visit in visits:
            # Get basic parameters
            viscode = visit['viscode']
            study = visit['study']
            filename = visit['filename']

            # Check viscode
            if args.viscode is not None and viscode != args.viscode:
                continue

            # Get image
            image = get_image(args, study, viscode, filename)

            # Get segmentation
            if args.method == 'graph':
                for structure in structures:
                    if args.structure_index is not None and adni.volume_names.index(structure) != args.structure_index:
                        continue

                    study_bl = visit['study_bl']
                    seg = get_segmentation_graphcut(study_bl, filename, structure)

                    print adni.INFO, 'Displaying {0} of subject {1}, ({2})'.format(structure, rid, viscode)
                    show_image(image, seg, lut=lut)
            else:
                if args.method == 'reg':
                    seg = get_segmentation_regbased(args, viscode, study, filename)
                elif args.method == 'long':
                    seg = get_segmentation_longitudinal(viscode, study, filename)
                elif args.method == 'cons':
                    seg = get_segmentation_consistent(sorted_viscodes.index(viscode), filename)

                print adni.INFO, 'Displaying segmentation of subject {0}, ({1})'.format(rid, viscode)
                show_image(image, seg, lut=lut)

                if args.method == 'cons':
                    remove_segmentation(seg)


def get_segmentation_regbased(args, viscode, study, filename):
    if viscode == 'bl':
        return adni.find_file(os.path.join(adni.data_folder, study, 'native/seg_138regions_baseline', 'EM-' + filename))
    else:
        return adni.find_file(os.path.join(adni.data_folder, study, 'native/seg_138regions_followup_' + args.trans + '_' + args.sx + 'mm', filename))


def get_segmentation_longitudinal(viscode, study, filename):
    if viscode == 'bl':
        return adni.find_file(os.path.join(adni.data_folder, study, 'native/seg_138regions_baseline', 'EM-' + filename))
    else:
        return adni.find_file(os.path.join('/vol/medic02/users/cl6311/data/ADNI/ADNI_followup/seg/EM/', 'EM-' + filename))


def get_segmentation_consistent(index, filename):
    print adni.INFO, 'Extracting slice {0} from {1}...'.format(index, filename)

    seg_4d = os.path.join('/vol/medic02/users/cl6311/ADNI4D/results', filename[:15] + '_4D.nii.gz')
    if not os.path.isfile(seg_4d):
        print adni.WARNING, 'Segmentation not found {0}'.format(seg_4d)
        return None
    else:
        seg_temp = os.path.join('/tmp/', filename)
        call([EXEC_EXTRACT, seg_4d, seg_temp, str(index)])

        return seg_temp


def get_segmentation_graphcut(study_bl, filename, structure):
    seg_dir = os.path.join(adni.data_folder, study_bl, 'native/seg_138regions_graphcut')
    postfix = '_' + structure.replace(' ', '_')
    return os.path.join(seg_dir, filename.replace('.nii.gz', postfix + '.nii.gz'))


def get_sorted_viscodes(cur, rid):
    import re
    import numpy as np

    cur.execute("SELECT viscode FROM adnimerge WHERE rid = " + str(rid))
    visits = cur.fetchall()
    scantimes = []
    viscodes = []
    for visit in visits:
        viscode = visit['viscode']
        if viscode == 'bl':
            scantime = 0
        elif re.match('m[0-9][0-9]', viscode):
            scantime = int(viscode[1:])
        else:
            print adni.ERROR, 'Invalid viscode', rid, viscode
            break

        scantimes.append(scantime)
        viscodes.append(viscode)

    scantimes = np.array(scantimes)
    viscodes = np.array(viscodes)

    indices = np.argsort(scantimes)
    return viscodes[indices].tolist()


def get_image(args, study, viscode, filename):
    if args.method in ('cons', 'graph') and viscode != 'bl':
        return os.path.join(adni.data_folder, study, 'baseline_linear/images_unstripped', filename)
    else:
        return os.path.join(adni.data_folder, study, 'native/images_unstripped', filename)


def show_image(image, segmentation, lut=None):
    if image is None or not os.path.isfile(image):
        print adni.ERROR, 'Image not found:', image
    if segmentation is None or not os.path.isfile(segmentation):
        print adni.ERROR, 'Segmentation not found:', segmentation
    else:
        print adni.INFO, 'Image:       ', image
        print adni.INFO, 'Segmentation:', segmentation
        if lut is None:
            call([EXEC_RVIEW, image, '-seg', segmentation, '-labels'])
        else:
            call([EXEC_RVIEW, image, '-seg', segmentation, '-lut', lut])


def remove_segmentation(segmentation):
    if segmentation is not None and segmentation[:5] == '/tmp/' and os.path.isfile(segmentation):
        call([EXEC_REMOVE, segmentation])

if __name__ == '__main__':
    main()
