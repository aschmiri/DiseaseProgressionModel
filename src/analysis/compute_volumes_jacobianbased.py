#! /usr/bin/env python
# print __doc__
import os.path
import argparse
import csv
import sqlite3
from subprocess import check_output
from src.common import adni_tools as adni


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--min_scans', dest='min_scans', type=int, default=6, help='the minimal number of scans per subject')
    parser.add_argument('-b', '--binary', dest='binary', action='store_true', help='select if binary segmentations should be used')
    a = parser.parse_args()

    exec_volumes = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/asr_compute_volume_probabilistic'
    exec_factors = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/asr_affineScaling'

    out_file = os.path.join(adni.project_folder, 'lists/volumes_jacobianbased.csv')

    # folder_lin_mni = adni.data_folder + '/ALL/MNI152_linear/dof'
    folder_lin_bl = adni.data_folder + '/ALL/baseline_linear/dof'
    folder_nonlin_bl = adni.data_folder + '/ALL/baseline_sym_10mm_after_linear/dof'

    #
    # Read csv file with metadata
    #
    con = sqlite3.connect(os.path.join(adni.project_folder, 'lists', 'adni.db'))
    con.row_factory = sqlite3.Row

    cur = con.cursor()
    cur.execute("SELECT rid FROM adnimerge WHERE viscode = 'bl'")
    rows = cur.fetchall()
    rids = [row['rid'] for row in rows]

    with open(out_file, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['RID', 'VISCODE', 'DX.scan'] + adni.volume_names)

        for rid in rids:
            cur.execute("SELECT iid, viscode, study, diagnosis, filename FROM adnimerge WHERE rid = " + str(rid))
            scans = cur.fetchall()
            if len(scans) >= a.min_scans:
                for scan in scans:
                    viscode = scan['viscode']
                    study = scan['study']
                    filename = scan['filename']
                    diagnosis = scan['diagnosis']

                    if viscode == 'bl':
                        bl_base = filename
                        bl_lin = None
                        bl_nonlin = None
                        bl_factor = 1
                    else:
                        cur.execute("SELECT filename FROM adnimerge WHERE rid = " + str(rid) + " AND viscode = 'bl'")
                        bl_base = cur.fetchall()[0]['filename']
                        bl_lin = os.path.join(folder_lin_bl.replace('ALL', study), filename.replace('.nii.gz', '.dof.gz'))
                        bl_nonlin = os.path.join(folder_nonlin_bl.replace('ALL', study), filename.replace('.nii.gz', '.dof.gz'))
                        if not os.path.isfile(bl_nonlin):
                            print 'ERROR: File not found:', bl_nonlin
                            break
                        if not os.path.isfile(bl_lin):
                            print 'ERROR: File not found:', bl_lin
                            break
                        bl_factor = float(check_output([exec_factors, bl_lin]))

                    # Get folder with the binary probablistic segmentations
                    if a.binary:
                        seg = adni.find_file(os.path.join(adni.data_folder, study, 'native/seg_138regions_baseline', 'EM-' + bl_base))
                        if seg is None or not os.path.isfile(seg):
                            print 'ERROR: File not found:', seg
                            break
                    else:
                        if study == 'ADNI1':
                            seg_folder_base = '/vol/medic02/users/cl6311/data/ADNI/ADNI1_baseline_1-5T/seg/prob_EM/'
                        elif study == 'ADNI2':
                            seg_folder_base = '/vol/medic02/users/cl6311/data/ADNI/ADNI2_baseline_3T/seg/prob_EM/'
                        elif study == 'ADNIGO':
                            seg_folder_base = '/vol/medic02/users/cl6311/data/ADNI/ADNIGO_baseline_3T/seg/prob_EM/'
                        else:
                            print 'ERROR: Invalid study'
                            break

                        seg = seg_folder_base + bl_base.replace('.nii.gz', '/')
                        if not os.path.isdir(seg):
                            print 'ERROR: Folder not found:', seg
                            break

                    # Get volumes of 138 objects
                    print '--------------------'
                    print 'bl ', bl_base
                    print viscode, filename
                    print 'seg', seg

                    if viscode == 'bl':
                        volumes = check_output([exec_volumes, seg])
                    else:
                        print 'def', bl_nonlin
                        volumes = check_output([exec_volumes, seg, bl_nonlin])
                    volumes = [float(vol) * bl_factor for vol in volumes.split(',')]
                    volumes.pop(0)
                    if len(volumes) != 138:
                        print 'ERROR:', len(volumes), 'volumes read for', seg
                    else:
                        writer.writerow([str(rid), viscode, diagnosis] + volumes)
                        csvfile.flush()
                        print rid, viscode, diagnosis, bl_factor, ':', volumes


if __name__ == '__main__':
    main()
