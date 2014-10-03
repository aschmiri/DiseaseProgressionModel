#! /bin/bash
folder="/vol/medic01/users/aschmidt/Dropbox/Documents/Einreichungen/NeuroImage 2014/figs/"
#folder="/Users/aschmiri/Dropbox/Documents/Einreichungen/NeuroImage 2014/figs/"

#
# Model generation
#
#biomarker="CDRSB"
#fitting/vgam_plot_curves.py cog -b "${biomarker}" --no_densities --no_extrapolation --no_sample_lines --no_model --output_file "${folder}model_generation_align.pdf"
#fitting/vgam_plot_curves.py cog -b "${biomarker}" --no_densities --no_extrapolation --no_sample_lines --output_file "${folder}model_generation_fit.pdf"
#fitting/vgam_plot_curves.py cog -b "${biomarker}" --no_densities --plot_mu --output_file "${folder}model_generation_extrapolate.pdf"
#fitting/vgam_plot_curves.py cog -b "${biomarker}" --only_densities --output_file "${folder}model_generation_densities.pdf"

#
# Synth models
#
#for b in hipp brain cdrsb mmse
#do
#  biomarker="synth_${b}"
#  output_file="${folder}synthetic_model_${b}.pdf"
#  fitting/vgam_plot_curves.py synth -b "${biomarker}" --plot_synth_model --no_densities --no_points --no_sample_lines --no_model --output_file "${output_file}"
#done

#
# Synth sampling
#
#biomarker="synth_cdrsb"
#for sampling in "uniform" "triangular" "longitudinal"
#do
#  output_file="${folder}synthetic_sampling_${sampling}.pdf"
#  fitting/vgam_generate_synth_data.py -b "${biomarker}" -n 100 --sampling ${sampling}
#	 fitting/vgam_estimate_curves.py synth -b "${biomarker}"
#  fitting/vgam_plot_curves.py synth -b "${biomarker}" --plot_synth_model --no_densities --no_sample_lines --no_extrapolation --no_model --points_alpha 1.0 --output_file "${output_file}"
#done

output_file="${folder}eval_synth_1.pdf"
#fitting/vgam_evaluate_synth_model.py -e ex1 -b synth_brain synth_hipp --experiment_range 100 1000 100 --number_of_runs 50 --output_file "${output_file}"
#fitting/vgam_evaluate_synth_model.py -e ex1 --output_file ${output_file}
output_file="${folder}eval_synth_2.pdf"
#fitting/vgam_evaluate_synth_model.py -e ex2 --output_file "${output_file}"
output_file="${folder}eval_synth_3.pdf"
#fitting/vgam_evaluate_synth_fitting_for_samplings.py --output_file "${output_file}"
output_file="${folder}eval_synth_4.pdf"
#fitting/vgam_evaluate_synth_fitting.py --output_file "${output_file}"

#
# Synth sampling
#
for method in cog mbl long all
do
  i=1
  for visits in "bl" "bl m12" "bl m12 m24"
  do
    output_file="${folder}eval_fitting_${method}_${i}.pdf"
    fitting/vgam_fit_subjects.py ${method} ${visits} --out "${output_file}"
    i=`expr $i + 1`   
  done
done


#
# Model examples
#
#for biomarker in "Right Hippocampus" "MMSE" "D1" "D2" "P_D1D2" 
#do
#  biomarker_str=${biomarker/ /_}
#  output_file="${folder}model_${biomarker_str}.pdf"
#  fitting/vgam_plot_curves.py all -b "${biomarker}" --output_file "${output_file}"
#done
