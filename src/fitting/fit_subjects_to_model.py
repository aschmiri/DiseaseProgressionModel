#! /usr/bin/env python2.7
import os.path
import csv
import re
import numpy as np
from common import adni_tools as adni


def main():
    load_model_data()


def exponential(t, r, r2, lower):
    y = lower + r * np.exp(r2 * t)
    return y


def load_model_data():
    global model_names
    global model_params

    model_file = os.path.join(adni.project_folder, 'lists/model_parameter.csv')
    with open(model_file, 'rb') as csvfile:
        rows = csv.DictReader(csvfile)
        for row in rows:
            # Get rid and age
            name = row['BIOMARKER']
            popt = np.array(float(row['PAR_R1']),
                            float(row['PAR_R2']),
                            float(row['PAR_LOWER']))

            model_names.append(name)
            model_params.append(popt)


def get_trajectory(biomarker_name, subject_rid):
    traj_x = []
    traj_y = []
    traj_d = []

    # Read data from CSV file
    data_file = os.path.join(adni.project_folder, 'lists/volumes_segbased_sym_5mm.csv')
    with open(data_file, 'rb') as csvfile:
        rows = csv.DictReader(csvfile)
        for row in rows:
            # Get rid and age
            rid = int(row['RID'])
            if rid == subject_rid:
                # Get scan time
                viscode = row['VISCODE']
                if viscode == 'bl':
                    scan_time = 0
                elif re.match('m[0-9][0-9]', viscode):
                    scan_time = int(viscode[1:])
                else:
                    print 'ERROR: Invalid viscode:', viscode
                    break

                # Get diagnosis
                dx_str = row['DX.scan']
                if dx_str == 'AD':
                    dx = 1.0
                elif dx_str == 'MCI':
                    dx = 0.5
                elif dx_str == 'CN':
                    dx = 0.0
                else:
                    print 'ERROR: Invalid diagnosis:', viscode
                    break

                # Get and normalise volumes
                value = adni.safe_cast(row[biomarker_name])
                if biomarker_name not in adni.cog_score_names:
                    value = value / float(row['FactorMNI'])

                if value is not None:
                    traj_x.append(scan_time)
                    traj_y.append(value)
                    traj_d.append(dx)

    # Sort data
    traj_x = np.array(traj_x)
    traj_y = np.array(traj_y)
    traj_d = np.array(traj_d)

    args = np.argsort(traj_x)
    traj_x = traj_x[args]
    traj_y = traj_y[args]
    traj_d = traj_d[args]

    return traj_x, traj_y, traj_d


if __name__ == '__main__':
    main()
