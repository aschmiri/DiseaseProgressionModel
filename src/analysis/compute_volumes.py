#! /usr/bin/env python2.7
import os.path
import csv
import sqlite3
from subprocess import check_output
from common import adni_tools as adni


def main():
    exec_volumes = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_compute_volume'
    exec_factors = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/affineScaling'

    data_file = os.path.join(adni.project_folder, 'lists/query_ADNI.csv')
    out_file = os.path.join(adni.project_folder, 'lists/volumes.csv')
    folder = adni.data_folder + '/ALL/native/images_unstripped'

    folder_lin_mni = adni.data_folder + '/ALL/MNI152_linear/dof'
    folder_lin_bl = adni.data_folder + '/ALL/baseline_linear/dof'

    #
    # Setup database
    #
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    cur.execute("CREATE TABLE summary (iid INTEGER PRIMARY KEY, \
                                        rid INTEGER, \
                                        factor REAL, \
                                        study TEXT, \
                                        viscode TEXT, \
                                        scandate TEXT, \
                                        diagnosis TEXT, \
                                        segname TEXT)")

    #
    # Read csv file with metadata
    #
    bl_base = None
    with open(data_file, 'rb') as csvfile:
        entries = csv.DictReader(csvfile)
        for entry in entries:
            study = entry['COLPROT']
            image = os.path.join(folder.replace('ALL', study), os.path.basename(entry['Files']))
            if os.path.isfile(image):
                image_base = os.path.basename(image)
                viscode = entry['VISCODE']
                fs = entry['MagStrength']

                if (study == 'ADNI1' and float(fs) < 2.0) or (study != 'ADNI1' and float(fs) > 2.0):
                    if viscode == 'bl':
                        seg = os.path.join(adni.data_folder, study, 'native/seg_138regions_baseline', 'EM-' + image_base)

                        bl_base = image_base
                        dof = os.path.join(folder_lin_mni.replace('ALL', study), bl_base.replace('.nii.gz', '.dof.gz'))
                        if os.path.isfile(dof):
                            factor = float(check_output([exec_factors, dof]))
                        else:
                            factor = -1
                    else:
                        seg = os.path.join(adni.data_folder, study, 'native/seg_138regions_followup_warped', image_base)

                        bl_dof = os.path.join(folder_lin_mni.replace('ALL', study), bl_base.replace('.nii.gz', '.dof.gz'))
                        fu_dof = os.path.join(folder_lin_bl.replace('ALL', study), image_base.replace('.nii.gz', '.dof.gz'))
                        if os.path.isfile(bl_dof) and os.path.isfile(fu_dof):
                            bl_factor = float(check_output([exec_factors, bl_dof]))
                            fu_factor = float(check_output([exec_factors, fu_dof]))
                            factor = bl_factor * fu_factor
                        else:
                            factor = -1

                    if os.path.isfile(seg) and factor > 0:
                        iid = entry['ImageUID']
                        rid = entry['RID']
                        dx = entry['DX.bl']
                        scandate = entry['ScanDate']

                        cur.execute("INSERT INTO summary VALUES (?,?,?,?,?,?,?,?)",
                                    (iid, rid, factor, study, viscode, scandate, dx, seg))

    #
    # Get database cursor
    #
    cur.execute("SELECT iid, rid, factor, study, viscode, scandate, diagnosis, segname FROM summary")
    rows = cur.fetchall()

    #
    # Write csv file with volumes
    #
    print 'Found', str(len(rows)), 'image pairs...'
    with open(out_file, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['RID', 'VISCODE', 'DX.bl', 'Factor'] + adni.volume_names)
        for row in rows:
            rid = row['rid']
            factor = row['factor']
            viscode = row['viscode']
            diagnosis = row['diagnosis']
            segname = row['segname']

            # Get volumes of 138 objects
            volumes = check_output([exec_volumes, segname])
            volumes = [float(vol) for vol in volumes.split(',')]
            volumes.pop(0)
            if len(volumes) != 138:
                print 'ERROR:', len(volumes), 'volumes read for', os.path.basename(segname)
            else:
                writer.writerow([str(rid), viscode, diagnosis, factor] + volumes)
            print rid, factor, viscode, diagnosis, ':', volumes


if __name__ == '__main__':
    main()
