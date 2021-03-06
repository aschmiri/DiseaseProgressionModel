#! /usr/bin/env python2.7
import os.path
import argparse
import math
import csv
import random
from common import log as log
from common.datahandler import SynthDataHandler
from common.synthmodel import SynthModel


EXEC_VOLUMES = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_compute_volume'
EXEC_EXTRACT = '/vol/biomedic/users/cl6311/irtk_svn_workspace/irtk/build/bin/cl_extract_frame'
EXEC_FACTORS = '/vol/medic01/users/aschmidt/development/build_helvellyn/irtk-as12312/bin/asr_affineScaling'
EXEC_REMOVE = 'rm'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--biomarkers', default=None, nargs='+', help='name of the biomarkers for which data is to be generated')
    parser.add_argument('-s', '--sampling', type=str, choices=['uniform', 'triangular', 'longitudinal'], default='triangular', help='the type of sampling')
    parser.add_argument('-n', '--number_of_samples', type=int, default=1000, help='the number of subjects in the synthetic sample pool')
    parser.add_argument('--samples_per_subject', type=int, default=6, help='the number of samples per subject')
    parser.add_argument('--days_between_samples', type=int, default=365, help='the temporal offset between to samples of one subject')
    parser.add_argument('--rate_sigma', type=float, default=0.0, help='the standard deviation of the gaussian noise applied on the progression rate')
    parser.add_argument('--conversion_sigma', type=float, default=0.0, help='the standard deviation of the gaussian noise applied on the point of conversion')
    args = parser.parse_args()

    print log.INFO, 'Generating synthetic data with {0} samples...'.format(args.number_of_samples)
    # Determine structures to segment
    data_handler = SynthDataHandler()
    biomarkers = data_handler.get_biomarker_names()

    # Estimate volumes
    out_file = os.path.join(data_handler.get_data_folder(), 'synth_data.csv')
    with open(out_file, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['RID', 'VISCODE', 'DX.scan', 'progress', 'rate', 'offset'] + biomarkers)

        if args.sampling == 'longitudinal':
            # Get the number of subjects from the total number of samples
            number_of_rids = int(math.ceil(float(args.number_of_samples) / float(args.samples_per_subject)))

            # Get samples of each subject
            for rid in range(number_of_rids):
                # Get days between two samples
                rate = 1.0 if args.rate_sigma == 0.0 else random.gauss(1.0, args.rate_sigma)
                sample_range = rate * args.days_between_samples * args.samples_per_subject - 1
                offset = 0.0 if args.conversion_sigma == 0.0 else random.gauss(0.0, args.conversion_sigma)
                base_progress_measured = SynthModel.get_random_progress(sampling=args.sampling, sample_range=sample_range)
                base_progress_real = base_progress_measured + offset

                # Get cdfs for all biomarkers
                cdfs = {}
                for biomarker in biomarkers:
                    base_sample = SynthModel.get_distributed_value(biomarker, base_progress_real)
                    cdfs.update({biomarker: SynthModel.get_cumulative_probability(biomarker, base_progress_real, base_sample)})

                # Get samples with same noise for all biomarkers
                for sample_index in range(args.samples_per_subject):
                    sample_progress_measured = base_progress_measured + args.days_between_samples * sample_index
                    sample_progress_real = int(round(base_progress_real + rate * args.days_between_samples * sample_index))

                    dx = 'AD' if sample_progress_real > 0 else 'MCI'
                    viscode = 'bl' if sample_index == 0 else 'm{0}'.format(sample_index * 12)
                    values = []
                    for biomarker in biomarkers:
                        values.append(SynthModel.get_distributed_value(biomarker, sample_progress_real, cdf=cdfs[biomarker]))
                    writer.writerow([rid, viscode, dx, sample_progress_measured, rate, offset] + values)
        else:
            for rid in range(args.number_of_samples):
                progress = SynthModel.get_random_progress(sampling=args.sampling)

                dx = 'AD' if progress > 0 else 'MCI'
                viscode = 'bl'
                values = []
                for biomarker in biomarkers:
                    values.append(SynthModel.get_distributed_value(biomarker, progress))
                writer.writerow([rid, viscode, dx, progress, '1.0', '0.0'] + values)


if __name__ == '__main__':
    main()
