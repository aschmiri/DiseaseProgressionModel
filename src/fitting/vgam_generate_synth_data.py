#! /usr/bin/env python2.7
import os.path
import argparse
import math
import csv
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import SynthDataHandler
from vgam.synthmodel import SynthModel


EXEC_VOLUMES = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_compute_volume'
EXEC_EXTRACT = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_extract_frame'
EXEC_FACTORS = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/asr_affineScaling'
EXEC_REMOVE = 'rm'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--biomarkers_name', default=None, nargs='+', help='name of the biomarkers for which data is to be generated')
    parser.add_argument('-s', '--sampling', type=str, choices=['uniform', 'triangular', 'longitudinal'], default='triangular', help='the type of sampling')
    parser.add_argument('-n', '--number_of_samples', type=int, default=400, help='the number of subjects in the synthetic sample pool')
    parser.add_argument('--samples_per_subject', type=int, default=6, help='the number of samples per subject')
    parser.add_argument('--temporal_offset', type=int, default=365, help='the temporal offset between to samples of one subject')
    args = parser.parse_args()

    print log.INFO, 'Generating synthetic data with {0} samples...'.format(args.number_of_samples)
    out_file = os.path.join(adni.project_folder, 'lists/synthdata.csv')

    # Determine structures to segment
    data_handler = SynthDataHandler(args)
    biomarkers = data_handler.get_biomarker_set()

    # Estimate volumes
    with open(out_file, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['RID', 'VISCODE', 'DX.scan', 'progress'] + biomarkers)

        if args.sampling == 'longitudinal':
            sample_range = args.samples_per_subject * args.temporal_offset - 1
            number_of_rids = int(math.ceil(float(args.number_of_samples) / float(args.samples_per_subject)))
            for rid in range(number_of_rids):
                base_progress = SynthModel.get_random_progress(sampling=args.sampling, sample_range=sample_range)

                # Get cdfs for all biomarkers
                cdfs = {}
                for biomarker in biomarkers:
                    base_sample = SynthModel.get_distributed_value(biomarker, base_progress)
                    cdfs.update({biomarker: SynthModel.get_cumulated_probability(biomarker, base_progress, base_sample)})

                # Get samples with same noise for all biomarkers
                for sample in range(args.samples_per_subject):
                    sample_progress = base_progress + sample * args.temporal_offset

                    dx = 'AD' if sample_progress > 0 else 'MCI'
                    viscode = 'bl' if sample == 0 else 'm{0}'.format(sample * 12)
                    values = []
                    for biomarker in biomarkers:
                        values.append(SynthModel.get_distributed_value(biomarker, sample_progress, cdf=cdfs[biomarker]))
                    writer.writerow([rid, viscode, dx, sample_progress] + values)
        else:
            for rid in range(args.number_of_samples):
                progress = SynthModel.get_random_progress(sampling=args.sampling)

                dx = 'AD' if progress > 0 else 'MCI'
                viscode = 'bl'
                values = []
                for biomarker in biomarkers:
                    values.append(SynthModel.get_distributed_value(biomarker, progress))
                writer.writerow([rid, viscode, dx, progress] + values)


if __name__ == '__main__':
    main()
