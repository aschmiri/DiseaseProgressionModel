#! /bin/bash
folder="/vol/medic01/users/aschmidt/Dropbox/Documents/Einreichungen/NeuroImage 2014/figs/"
#folder="/Users/aschmiri/Dropbox/Documents/Einreichungen/NeuroImage 2014/figs/"

#
# Model generation
#
#biomarker="CDRSB"
#fitting/vgam_plot_curves.py -b "${biomarker}" --no_densities --no_extrapolation --no_sample_lines --no_model --plot_file "${folder}model_generation_align.pdf"
#fitting/vgam_plot_curves.py -b "${biomarker}" --no_densities --no_extrapolation --no_sample_lines --plot_file "${folder}model_generation_fit.pdf"
#fitting/vgam_plot_curves.py -b "${biomarker}" --no_densities --plot_mu --plot_file "${folder}model_generation_extrapolate.pdf"
#fitting/vgam_plot_curves.py -b "${biomarker}" --only_densities --plot_file "${folder}model_generation_densities.pdf"

#
# Synth models
#
#for b in hipp cdrsb mmse
#do
#  biomarker="synth_${b}"
#  output_file="${folder}synthetic_model_${b}.pdf"
#  fitting/vgam_generate_synth_data.py -b "${biomarker}" -n 1000 --sampling longitudinal
#  fitting/vgam_estimate_curves.py -m synth -b "${biomarker}"
#  fitting/vgam_plot_curves.py -m synth -b "${biomarker}" --plot_synth_model --no_densities --no_points --no_sample_lines --plot_file "${output_file}"
#done

#
# Synth sampling
#
#biomarker="synth_cdrsb"
#for sampling in "longitudinal" "triangular" "uniform"
#do
#  output_file="${folder}synthetic_sampling_${sampling}.pdf"
#  fitting/vgam_generate_synth_data.py -b "${biomarker}" -n 100 --sampling ${sampling}
#  fitting/vgam_estimate_curves.py -m synth -b "${biomarker}" --recompute_models
#  fitting/vgam_plot_curves.py -m synth -b "${biomarker}" --plot_synth_model --no_densities --no_sample_lines --no_extrapolation --no_model --points_alpha 0.7 --plot_file "${output_file}"
#done

#
# Synth evaluation
#
#output_file="${folder}eval_synth_1.pdf"
#fitting/vgam_evaluate_synth_model.py -e ex1 -b synth_cdrsb --sample_numbers_range 100 2000 100 --number_of_runs 100 --plot_file "${output_file}"
#output_file="${folder}eval_synth_2.pdf"
#fitting/vgam_evaluate_synth_model.py -e ex2 --plot_file "${output_file}"
#output_file="${folder}eval_synth_sigma_rates.pdf"
#fitting/vgam_evaluate_synth_model.py -e ex3 --plot_file "${output_file}"
#output_file="${folder}eval_synth_sigma_conv.pdf"
#fitting/vgam_evaluate_synth_model.py -e ex4 --plot_file "${output_file}"

#output_file="${folder}eval_synth_3.pdf"
#fitting/vgam_evaluate_synth_fitting.py -e ex1 --plot_file "${output_file}"
output_file="${folder}eval_synth_4.pdf"
fitting/vgam_evaluate_synth_fitting.py -e ex2 --plot_file "${output_file}"


#
# Model examples
#
#for biomarker in "Right Hippocampus" "MMSE" "D1" "D2" "P_D1D2" 
#do
#  biomarker_str=${biomarker/ /_}
#  output_file="${folder}model_${biomarker_str}.pdf"
#  fitting/vgam_plot_curves.py -b "${biomarker}" --plot_file "${output_file}"
#done

#
# Progression estimates
#
#for method in cog mbl long img all
#do
#  i=1
#  for visits in "bl" "bl m12" "bl m12 m24"
#  do
#    output_file="${folder}eval_fitting_${method}_${i}.pdf"
#    fitting/vgam_estimate_progress.py ${visits} -m ${method} --plot_file "${output_file}"
#    i=`expr $i + 1`   
#  done
#done

#
# Biomarker predictions
#
#for biomarker in "MMSE" "CDRSB" "FAQ" "Right Hippocampus"
#do
#  biomarker_str=${biomarker/ /_}
#  output_file="${folder}predict_${biomarker_str}.pdf"
#  fitting/vgam_compare_biomarker_predictions.py  -p "${biomarker}" --plot_file "${output_file}"
#done

