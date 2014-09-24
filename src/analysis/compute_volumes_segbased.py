#! /usr/bin/env python2.7
import os.path
import argparse
import csv
import sqlite3
from subprocess import call
from subprocess import check_output
from common import log as log
from common import adni_tools as adni


EXEC_VOLUMES = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_compute_volume'
EXEC_EXTRACT = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_extract_frame'
EXEC_FACTORS = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/asr_affineScaling'
EXEC_REMOVE = 'rm'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('method', choices=['reg', 'long', 'cons', 'graph', 'meta'])
    parser.add_argument('-m', '--min_scans', type=int, default=0, help='the minimal number of scans per subject')
    parser.add_argument('--trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic (regbased only)')
    parser.add_argument('--spacing', type=str, default='5', help='the transformation spacing (regbased only)')
    parser.add_argument('--all_structures', action='store_true', default=False, help='segment all 138 structures (graphcuts only)')
    args = parser.parse_args()

    # Set folders for factor computation
    global folder_lin_mni
    global folder_lin_bl
    folder_lin_mni = adni.data_folder + '/ALL/MNI152_linear/dof'
    folder_lin_bl = adni.data_folder + '/ALL/baseline_linear/dof'

    out_file_ending = {'reg': 'lists/volumes_segbased_' + args.trans + '_' + args.spacing + 'mm.csv',
                       'long': 'lists/volumes_segbased_longitudinal.csv',
                       'cons': 'lists/volumes_segbased_consistent.csv',
                       'graph': 'lists/volumes_segbased_graphcut.csv',
                       'meta': 'lists/metadata.csv'}
    out_file = os.path.join(adni.project_folder, out_file_ending[args.method])

    # Determine structures to segment
    if args.method == 'meta':
        structures = []
    elif args.method == 'graph' and not args.all_structures:
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

    # Estimate volumes
    with open(out_file, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['RID', 'VISCODE', 'DX.scan', 'AGE.scan', 'ScanDate'] +
                        adni.cog_score_names +
                        ['FactorMNI', 'FactorBL'] +
                        structures)

        for rid in rids:
            # Get MNI scaling for RID
            mni_factor = get_mni_factor(cur, rid)

            # Get sorted viscodes (only needed for consistent 4D segmentation)
            sorted_viscodes = get_sorted_viscodes(cur, rid)

            # Get all visits for RID
            cur.execute("SELECT iid, viscode, study, study_bl, diagnosis, filename, age, scandate \
                          FROM adnimerge WHERE rid = " + str(rid))
            visits = cur.fetchall()
            if len(visits) >= args.min_scans:
                for visit in visits:
                    # Get basic parameters
                    viscode = visit['viscode']
                    study = visit['study']
                    study_bl = visit['study_bl']
                    filename = visit['filename']
                    diagnosis = visit['diagnosis']
                    age = visit['age']
                    scandate = visit['scandate']

                    # Get cognitive scores
                    bl_factor = get_bl_factor(cur, rid, viscode, study, filename)

                    # Get cognitive scores
                    cog_scores = get_cog_scores(cur, rid, viscode)

                    # Get volumes
                    if args.method == 'reg':
                        volumes = get_volumes_regbased(args, rid, viscode, study, filename)
                    elif args.method == 'long':
                        volumes = get_volumes_longitudinal(rid, viscode, study, filename)
                    elif args.method == 'cons':
                        volumes = get_volumes_consistent(rid, sorted_viscodes.index(viscode), study, filename)
                    elif args.method == 'graph':
                        volumes = get_volumes_graphcut(rid, viscode, study_bl, filename, structures)
                    else:
                        volumes = []

                    if len(volumes) != len(structures):
                        print log.WARNING, '{0} volumes read for subject {1} ({2})'.format(len(volumes), rid, viscode)
                    else:
                        writer.writerow([str(rid), viscode, diagnosis, age, scandate] +
                                        cog_scores +
                                        [mni_factor, bl_factor] +
                                        volumes)
                        csvfile.flush()


def get_mni_factor(cur, rid):
    print log.INFO, 'Computing MNI scaling factor for subject {0}...'.format(rid)
    cur.execute("SELECT study, filename \
                FROM adnimerge WHERE rid = " + str(rid) + " AND viscode = 'bl'")
    bl_data = cur.fetchall()
    if len(bl_data) != 1:
        print log.WARNING, 'bl_data has wrong size ({0}), using first entry.'.format(len(bl_data))

    study = bl_data[0]['study']
    filename = bl_data[0]['filename']

    mni_dof = os.path.join(folder_lin_mni.replace('ALL', study), filename.replace('.nii.gz', '.dof.gz'))
    if not os.path.isfile(mni_dof):
        print log.ERROR, 'File not found:', mni_dof
        return None
    else:
        return float(check_output([EXEC_FACTORS, mni_dof]))


def get_bl_factor(cur, rid, viscode, study, filename):
    print log.INFO, 'Computing scaling factor for subject {0} ({1})...'.format(rid, viscode)
    if viscode == 'bl':
        return 1.0
    else:
        bl_dof = os.path.join(folder_lin_bl.replace('ALL', study), filename.replace('.nii.gz', '.dof.gz'))
        if not os.path.isfile(bl_dof):
            print log.ERROR, 'File not found:', bl_dof
            return None
        else:
            return float(check_output([EXEC_FACTORS, bl_dof]))


def get_cog_scores(cur, rid, viscode):
    # TODO: make consistent with adni.cogscore_names
    print log.INFO, 'Querying cognitive scores for subject {0} ({1})...'.format(rid, viscode)
    cur.execute("SELECT mmse, cdrsb, adas11, adas13, faq \
                FROM adnimerge JOIN adnimerge_cog USING (rid, viscode) \
                WHERE rid = " + str(rid) + " AND viscode = '" + viscode + "'")
    cog_data = cur.fetchall()
    if len(cog_data) != 1:
        print log.ERROR, 'cog_data has wrong size:', len(cog_data)

    mmse = adni.safe_cast(cog_data[0]['mmse'])
    cdrsb = adni.safe_cast(cog_data[0]['cdrsb'])
    adas11 = adni.safe_cast(cog_data[0]['adas11'])
    adas13 = adni.safe_cast(cog_data[0]['adas13'])
    faq = adni.safe_cast(cog_data[0]['faq'])

    return [mmse, cdrsb, adas11, adas13, faq]


def get_volumes_regbased(args, rid, viscode, study, filename):
    print log.INFO, 'Reading volumes for subject {0} ({1})...'.format(rid, viscode)
    if viscode == 'bl':
        seg = adni.find_alternative_file(os.path.join(adni.data_folder, study, 'native/seg_138regions_baseline', 'EM-' + filename))
    else:
        seg = adni.find_alternative_file(os.path.join(adni.data_folder, study, 'native/seg_138regions_followup_' + args.trans + '_' + args.spacing + 'mm', filename))

    # Get volumes of the cortical structures
    if seg is None:
        return []
    else:
        volumes = check_output([EXEC_VOLUMES, seg])
        volumes = [float(vol) for vol in volumes.split(',')]
        volumes.pop(0)

        return volumes


def get_volumes_longitudinal(rid, viscode, study, filename):
    print log.INFO, 'Reading volumes for subject {0} ({1})...'.format(rid, viscode)
    #if viscode == 'bl':
    #    seg = adni.find_alternative_file(os.path.join(adni.data_folder, study, 'native/seg_138regions_baseline', 'EM-' + filename))
    #else:
    #    seg = adni.find_alternative_file(os.path.join('/vol/medic02/users/cl6311/data/ADNI/ADNI_followup/seg/EM/', 'EM-' + filename))
    seg = adni.find_alternative_file(os.path.join('/vol/medic02/users/cl6311/ADNI_condor/results/','MALPEM-' + filename))
    
    # Get volumes of the cortical structures
    if seg is None:
        return []
    else:
        volumes = check_output([EXEC_VOLUMES, seg])
        volumes = [float(vol) for vol in volumes.split(',')]
        volumes.pop(0)

        return volumes


def get_volumes_consistent(rid, index, study, filename):
    print log.INFO, 'Reading volumes for subject {0} (index {1})...'.format(rid, index)

    seg_4d = os.path.join('/vol/medic02/users/cl6311/ADNI4D/results', filename[:15] + '_4D.nii.gz')
    if not os.path.isfile(seg_4d):
        print log.WARNING, 'Segmentation not found {0}'.format(seg_4d)
        return []
    else:
        seg_temp = os.path.join('/tmp/', filename)
        call([EXEC_EXTRACT, seg_4d, seg_temp, str(index)])

        if not os.path.isfile(seg_temp):
            print log.WARNING, 'Error occurred while extracting slice {0} from {1}'.format(slice, seg_4d)
            return []
        else:
            volumes = check_output([EXEC_VOLUMES, seg_temp])
            volumes = [float(vol) for vol in volumes.split(',')]
            volumes.pop(0)

            if os.path.isfile(seg_temp):
                call([EXEC_REMOVE, seg_temp])

            return volumes


def get_volumes_graphcut(rid, viscode, study_bl, filename, structures):
    # Get data for volume computation
    print log.INFO, 'Collecting volumes for subject {0} ({1})...'.format(rid, viscode)
    seg_dir = os.path.join(adni.data_folder, study_bl, 'native/seg_138regions_graphcut')
    volumes = []
    for structure in structures:
        postfix = '_' + structure.replace(' ', '_')

        seg = os.path.join(seg_dir, filename.replace('.nii.gz', postfix + '.nii.gz'))
        # Get volumes of 138 objects
        if not os.path.isfile(seg):
            print log.WARNING, 'Segmentation not found: {0}'.format(os.path.basename(seg))
            volumes.append(-1.0)
        else:
            output = check_output([EXEC_VOLUMES, seg])
            output_split = output.split(',')
            if len(output_split) == 1:
                print log.WARNING, 'No {0} found in {1}'.format(structure, os.path.basename(seg))
                volumes.append(0.0)
            elif len(output_split) == 2:
                volumes.append(float(output_split[1]))
            else:
                print log.WARNING, 'Corrupted output from {0} ({1})'.format(os.path.basename(seg), output_split)
                volumes.append(-1.0)

    return volumes


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
            print log.ERROR, 'Invalid viscode', rid, viscode
            break

        scantimes.append(scantime)
        viscodes.append(viscode)

    scantimes = np.array(scantimes)
    viscodes = np.array(viscodes)

    indices = np.argsort(scantimes)
    return viscodes[indices].tolist()


if __name__ == '__main__':
    main()
