#! /usr/bin/env python2.7
import os.path
import argparse
import csv
from common import log as log
from common import adni_tools as adni
from vgam.synthmodel import SynthModel


EXEC_VOLUMES = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_compute_volume'
EXEC_EXTRACT = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_extract_frame'
EXEC_FACTORS = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/asr_affineScaling'
EXEC_REMOVE = 'rm'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number_of_subjects', type=int, default=400, help='the number of subjects in the synthetic sample pool')
    parser.add_argument('-u', '--uniform_progression', action='store_true', help='use an uniform progression distribution')
    args = parser.parse_args()

    print log.INFO, 'Generating synthetic data with {0} samples...'.format(args.number_of_subjects)
    out_file = os.path.join(adni.project_folder, 'lists/synthdata.csv')

    # Determine structures to segment
    biomarkers = SynthModel.get_biomarker_names()

    # Estimate volumes
    with open(out_file, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['RID', 'VISCODE', 'DX.scan', 'progress'] + biomarkers)
        for rid in range(args.number_of_subjects):
            viscode = 'bl'
            progress = SynthModel.get_random_progress(uniform_progression=args.uniform_progression)

            dx = 'AD' if progress > 0 else 'MCI'
            values = []
            for biomarker in biomarkers:
                values.append(SynthModel.get_distributed_value(biomarker, progress))
            writer.writerow([rid, viscode, dx, progress] + values)


if __name__ == '__main__':
    main()
