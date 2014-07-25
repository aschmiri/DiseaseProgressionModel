#! /usr/bin/env python
# print __doc__
import argparse
import os.path
import joblib as jl
from subprocess import call
from src.common import adni_tools as adni


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO')
    parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    parser.add_argument('trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('-s', '--spacing', dest='sx', type=str, default='10')
    parser.add_argument('-n', '--nr_threads', dest='nr_threads', type=int, default=1)
    parser.add_argument('-f', '--forward', action='store_true', default=False, help='use forward registration from bl to fu')
    a = parser.parse_args()

    global dof_folder_lin
    global dof_folder_nonlin
    global dof_folder_followup
    data_folder = os.path.join(adni.data_folder, a.study)
    dof_folder_lin = os.path.join(data_folder, 'baseline_linear', 'dof')
    dof_folder_nonlin = os.path.join(data_folder, 'baseline_' + a.trans + '_' + a.sx + 'mm_after_linear', 'dof')
    dof_folder_followup = os.path.join(data_folder, 'followup_' + a.trans + '_' + a.sx + 'mm_after_linear', 'dof')

    global seg_folder_in
    global seg_folder_out
    seg_folder = os.path.join('/vol/medic01/users/aschmidt/projects/Data/ADNI/data', a.study, 'native')
    seg_folder_in = os.path.join(seg_folder, 'seg_138regions_baseline')
    seg_folder_out = adni.make_dir(seg_folder, 'seg_138regions_followup_' + a.trans + '_' + a.sx + 'mm')

    global baseline_files
    global followup_files
    baseline_folder = os.path.join(data_folder, 'native/images_unstripped')
    followup_folder = os.path.join(data_folder, 'native/images_unstripped')
    baseline_files, followup_files = adni.get_baseline_and_followup(baseline_folder, followup_folder, a.study, a.viscode)

    print 'Found', len(baseline_files), 'image pairs...'
    if a.forward:
        jl.Parallel(n_jobs=a.nr_threads)(jl.delayed(run_forward)(i, a.study) for i in range(len(baseline_files)))
    else:
        jl.Parallel(n_jobs=a.nr_threads)(jl.delayed(run_backward)(i, a.study) for i in range(len(baseline_files)))


def run_backward(index, study):
    target = baseline_files[index]
    target_base = os.path.basename(target)
    study_bl = adni.detect_study(target)
    source = followup_files[index]
    source_base = os.path.basename(source)

    dof_lin = os.path.join(dof_folder_lin, source_base.replace('.nii.gz', '.dof.gz'))
    dof_nonlin = os.path.join(dof_folder_nonlin, source_base.replace('.nii.gz', '.dof.gz'))

    target_seg = os.path.join(seg_folder_in.replace(study, study_bl), 'EM-' + target_base)
    target_seg = adni.find_file(target_seg)
    out_seg = os.path.join(seg_folder_out, source_base)

    if target_seg is not None:
        if not os.path.isfile(dof_nonlin):
            print 'DOF file', dof_nonlin, 'does not exists!'
        elif os.path.isfile(out_seg):
            print 'File', out_seg, 'already exists!'
        else:
            #
            # Run transform masks
            #
            print '--------------------'
            print 'Starting: transformation'
            print 'Target:    ', target
            print 'Source:    ', source
            print 'DOF lin:   ', dof_lin
            print 'DOF nonlin:', dof_nonlin
            print 'Seg in:    ', target_seg
            print 'Seg out:   ', out_seg

            dof_lin_ffd = out_seg.replace('.nii.gz', '_lin.dof.gz')
            dof_combined = out_seg.replace('.nii.gz', '_combined.dof.gz')

            call(['ffdcreate', dof_lin_ffd, '-dofin', dof_lin])
            call(['ffdcompose', dof_nonlin, dof_lin_ffd, dof_combined])
            call(['transformation', target_seg, out_seg, '-dofin', dof_combined, '-target', target, '-nn', '-matchInputType', '-invert'])

            call(['rm', dof_lin_ffd])
            call(['rm', dof_combined])

            print '--------------------'


def run_forward(index, study):
    target = followup_files[index]
    target_base = os.path.basename(target)
    source = baseline_files[index]
    source_base = os.path.basename(source)
    study_bl = adni.detect_study(source)

    dof = os.path.join(dof_folder_followup, 'baseline_to_' + target_base.replace('.nii.gz', '.dof.gz'))

    source_seg = os.path.join(seg_folder_in.replace(study, study_bl), 'EM-' + source_base)
    source_seg = adni.find_file(source_seg)
    out_seg = os.path.join(seg_folder_out, target_base)

    if source_seg is not None:
        if not os.path.isfile(dof):
            print 'DOF file', dof, 'does not exists!'
        elif os.path.isfile(out_seg):
            print 'File', out_seg, 'already exists!'
        else:
            #
            # Run transform masks
            #
            print '--------------------'
            print 'Starting: transformation'
            print 'Target:    ', target
            print 'Source:    ', source
            print 'DOF:       ', dof
            print 'Seg in:    ', source_seg
            print 'Seg out:   ', out_seg

            call(['transformation', source_seg, out_seg, '-dofin', dof, '-target', target, '-nn', '-matchInputType'])

            print '--------------------'
