#! /usr/bin/env python2.7
import os.path
import argparse
import csv
import numpy as np
from subprocess import call
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from common import log as log
from common import adni_tools as adni
from vgam.datahandler import SynthDataHandler
from vgam.synthmodel import SynthModel
from vgam.progressionmodel import ProgressionModel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--recompute_errors', action='store_true', help='recompute the errors of the models')
    parser.add_argument('--range', type=int, nargs=3, default=[100, 2000, 100], help='the range for the number of samples tested')
    parser.add_argument('--number_of_runs', type=int, default=5, help='the number of repeated runs')
    parser.add_argument('--number_of_progression_steps', type=int, default=10, help='the number of progression steps')
    parser.add_argument('--number_of_value_steps', type=int, default=1000, help='the number of value steps')
    parser.add_argument('--progression_range', type=int, default=2000, help='the width of progression range window')
    args = parser.parse_args()

    data_handler = SynthDataHandler()

    if args.recompute_errors:
        errors = compute_errors(args, data_handler)
        save_errors(args, data_handler, errors)
    else:
        errors = read_errors(args)

    plot_errors(args, data_handler, errors)
    analyse_errors(args, data_handler, errors)


def read_errors(args):
    evaluation_file = os.path.join(adni.eval_folder, 'eval_synth_model.csv')

    m = mlab.csv2rec(evaluation_file)
    biomarkers = m['biomarker']
    uniforms = m['uniform']
    experiments = np.sort([int(e) for e in m.dtype.names[2:]])

    errors = {}
    for experiment in experiments:
        errors.update({experiment: {}})
        measures = m[str(experiment)]
        for i in xrange(len(biomarkers)):
            biomarker = biomarkers[i]
            uniform = uniforms[i]
            if biomarker not in errors[experiment]:
                errors[experiment].update({biomarker: {}})
            errors[experiment][biomarker].update({uniform: measures[i]})

    return errors


def compute_errors(args, data_handler):
    errors = {}
    for num_samples in xrange(args.range[0], args.range[1], args.range[2]):
        errors.update({num_samples: evaluate_sample_number(args, data_handler, num_samples)})
    return errors


def evaluate_sample_number(args, data_handler, num_samples):
    print log.INFO, 'Evaluating synthetic models with {0} samples...'.format(num_samples)

    mean_errors = {}
    for biomarker in data_handler.get_biomarker_set():
        mean_errors.update({biomarker: {True: 0.0, False: 0.0}})

    for uniform in [True, False]:
        for _ in xrange(args.number_of_runs):
            generate_models(num_samples, uniform_progression=uniform)
            for biomarker in data_handler.get_biomarker_set():
                mean_errors[biomarker][uniform] += evaluate_model(args, data_handler, biomarker)

        for biomarker in data_handler.get_biomarker_set():
            mean_errors[biomarker][uniform] /= args.number_of_runs

    return mean_errors


def generate_models(num_samples, uniform_progression=False):
    folder = os.path.join(adni.project_folder, 'src', 'fitting')
    up_arg = '--uniform_progression' if uniform_progression else ''
    call('{0}/vgam_generate_synth_data.py {1} -n {2}'.format(folder, up_arg, num_samples), shell=True)
    call('{0}/vgam_estimate_curves.py synth'.format(folder), shell=True)


def evaluate_model(args, data_handler, biomarker):
    model_file = data_handler.get_model_file(biomarker)

    # Define progression steps
    pm = ProgressionModel(biomarker, model_file)
    progressions = np.linspace(-args.progression_range, args.progression_range, args.number_of_progression_steps)

    # Define value steps
    min_val = float('inf')
    max_val = float('-inf')
    for std in [-2.0, 2.0]:
        curve = pm.get_quantile_curve(progressions, std)
        min_val = min(min_val, np.min(curve))
        max_val = max(max_val, np.max(curve))
    values = np.linspace(min_val, max_val, args.number_of_value_steps)

    # Get mean error
    error = 0
    for progr in progressions:
        probs_fit = pm.get_density_distribution(values, progr)
        probs_model = [SynthModel.get_probability(biomarker, progr, v) for v in values]
        error += np.sum(np.abs(np.array(probs_fit) - np.array(probs_model)))
    error *= (values[1] - values[0]) / len(progressions)

    return error


def save_errors(args, data_handler, errors):
    evaluation_file = os.path.join(adni.eval_folder, 'eval_synth_model.csv')
    experiments = np.sort(errors.keys()).tolist()

    writer = csv.writer(open(evaluation_file, 'wb'), delimiter=',')
    writer.writerow(['biomarker', 'uniform'] + experiments)
    for biomarker in data_handler.get_biomarker_set():
        for uniform in [True, False]:
            error_curve = []
            for experiment in experiments:
                error_curve.append(errors[experiment][biomarker][uniform])
            writer.writerow([biomarker, uniform] + error_curve)


def plot_errors(args, data_handler, errors):
    plt.figure()
    linestyle = {True: '--', False: '-'}
    color = {'synth1': 'b', 'synth2': 'g', 'synth3': 'r', 'synth4': 'k'}

    experiments = np.sort(errors.keys())
    for biomarker in data_handler.get_biomarker_set():
        for uniform in [True, False]:
            error_curve = []
            for experiment in experiments:
                error_curve.append(errors[experiment][biomarker][uniform])
            plt.plot(experiments, error_curve,
                     linestyle=linestyle[uniform], color=color[biomarker],
                     label='{0}{1}'.format(biomarker, ' (unif)' if uniform else ''))

    plt.legend()
    plt.show()


def analyse_errors(args, data_handler, errors):
    experiments = np.sort(errors.keys())

    # Compute mean error uniform vs. triangular
    mean_difference = 0.0
    for biomarker in data_handler.get_biomarker_set():
        error_curve_uniform = []
        error_curve_triangu = []
        for experiment in experiments:
            error_curve_uniform.append(errors[experiment][biomarker][True])
            error_curve_triangu.append(errors[experiment][biomarker][False])

        mean_difference += np.mean(np.array(error_curve_triangu) -
                                   np.array(error_curve_uniform))
    mean_difference /= len(data_handler.get_biomarker_set())

    print log.RESULT, 'Mean difference uniform vs. triangular: {0}'.format(mean_difference)


if __name__ == '__main__':
    main()
