#! /usr/bin/env python2.7
import os.path
import argparse
from subprocess import call
from bmia.common import adni_tools as adni

EXEC_MODEL = 'model_generation'


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    parser.add_argument('diagnosis', type=str, help='the diagnosis, e.g. AD, MCI, CN, ...')
    parser.add_argument('-i', '--iteration', dest='iteration', type=int, default=0, help='the current iteration')
    parser.add_argument('-t', '--trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('--state_min', type=float, default=0, help='minimal state')
    parser.add_argument('--state_max', type=float, default=15, help='maximal state')
    parser.add_argument('--state_stepsize', type=float, default=0.25, help='step size for model generation')
    parser.add_argument('--postfix', type=str, default='')
    a = parser.parse_args()

    atlas_folder = os.path.join(adni.project_folder, 'atlas/model_' + str(a.iteration) + a.postfix)
    model_prefix = os.path.join(atlas_folder, 'model_' + a.trans + '_' + a.viscode + '_' + a.diagnosis + '_p')

    linear = True if a.iteration == 0 else False
    if linear:
        velo = os.path.join(atlas_folder, 'velo_' + a.trans + '_' + a.viscode + '_' + a.diagnosis + '.dof.gz')
        call([EXEC_MODEL, adni.mni_atlas, velo,
              '-min', str(a.state_min), '-max', str(a.state_max), '-stepsize', str(a.state_stepsize),
              '-saveimg', '-prefix', model_prefix, '-linear'])
    else:
        velo_file = os.path.join(atlas_folder, 'velo_' + a.trans + '_' + a.viscode + '_' + a.diagnosis + '.txt')
        call([EXEC_MODEL, adni.mni_atlas, velo_file,
              '-min', str(a.state_min), '-max', str(a.state_max), '-stepsize', str(a.state_stepsize),
              '-saveimg', '-prefix', model_prefix])


if __name__ == '__main__':
    main()
