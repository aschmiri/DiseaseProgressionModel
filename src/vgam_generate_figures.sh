#! /bin/bash
experiment=${1-"0"}

################################ 
###          MedIA           ###
################################
folder="/vol/medic01/users/aschmidt/Dropbox/Documents/Einreichungen/MedIA 2014/figs/"
#folder="/Users/aschmiri/Dropbox/Documents/Einreichungen/MedIA 2014/figs/"

#
# Model generation
#
if [ $experiment -eq "1" ] ; then
  biomarker="CDRSB"
  training/train_models.py
  training/plot_models.py -b "${biomarker}" --no_densities --no_extrapolation --no_sample_lines --no_model --plot_file "${folder}model_generation_align.pdf"
  training/plot_models.py -b "${biomarker}" --no_densities --no_extrapolation --no_sample_lines --plot_file "${folder}model_generation_fit.pdf"
  training/plot_models.py -b "${biomarker}" --no_densities --plot_file "${folder}model_generation_extrapolate.pdf"
  training/plot_models.py -b "${biomarker}" --only_densities --plot_file "${folder}model_generation_densities.pdf"
fi

#
# Synth models
#
if [ $experiment -eq "2" ] ; then
  for b in hipp cdrsb mmse
  do
    biomarker="synth_${b}"
    output_file="${folder}synthetic_model_${b}.pdf"
    synth/generate_data.py -b "${biomarker}" -n 1000 --sampling longitudinal
    training/train_models.py -m synth -b "${biomarker}"
    training/plot_models.py -m synth -b "${biomarker}" --plot_synth_model --no_densities --no_points --no_sample_lines --plot_file "${output_file}"
  done
fi

#
# Synth sampling
#
if [ $experiment -eq "3" ] ; then
  biomarker="synth_cdrsb"
  for sampling in "longitudinal" "triangular" "uniform"
  do
    output_file="${folder}synthetic_sampling_${sampling}.pdf"
    synth/generate_data.py -b "${biomarker}" -n 100 --sampling ${sampling}
    training/train_models.py -m synth -b "${biomarker}" --recompute_models
    training/plot_models.py -m synth -b "${biomarker}" --plot_synth_model --no_densities --no_sample_lines --no_extrapolation --no_model --points_alpha 0.7 --plot_file "${output_file}"
  done
fi

#
# Synth evaluation
#
if [ $experiment -eq "4" ] ; then
  output_file="${folder}eval_synth_num_samples.pdf"
  synth/evaluate_reconstruction.py -e ex1 -b synth_cdrsb --plot_file "${output_file}"
  output_file="${folder}eval_synth_sampling_strategy.pdf"
  synth/evaluate_reconstruction.py -e ex2 --plot_file "${output_file}"
  output_file="${folder}eval_synth_sigma_rates.pdf"
  synth/evaluate_reconstruction.py -e ex3 --plot_file "${output_file}"
  output_file="${folder}eval_synth_sigma_conv.pdf"
  synth/evaluate_reconstruction.py -e ex4 --plot_file "${output_file}"
fi

if [ $experiment -eq "5" ] ; then
  output_file="${folder}eval_synth_estimation.pdf"
  synth/evaluate_estimation.py -e ex1 --plot_file "${output_file}"
  output_file="${folder}eval_synth_available_data.pdf"
  synth/evaluate_estimation.py -e ex2 --plot_file "${output_file}"
fi


#
# Model examples
#
if [ $experiment -eq "6" ] ; then
  for biomarker in "Right Hippocampus" "MMSE" "D1" "D2"
  do
    biomarker_str=${biomarker/ /_}
    output_file="${folder}model_${biomarker_str}.pdf"
    training/plot_models.py -b "${biomarker}" --plot_file "${output_file}"
  done
fi

#
# Progression estimates
#
if [ $experiment -eq "7" ] ; then
  for method in  cog ml vol img all
  do
    i=1
    for visits in "bl" "bl m12" "bl m12 m24"
    do
      output_file="${folder}eval_fitting_${method}_${i}.pdf"
      estimation/estimate_progressions.py ${visits} -m ${method} --consistent_data --no_analysis --plot_file "${output_file}"
      i=`expr $i + 1`   
    done
  done
fi

#
# Biomarker predictions
#
if [ $experiment -eq "8" ] ; then
  for biomarker in "MMSE" "CDRSB" "FAQ" "Right Hippocampus"
  do
    biomarker_str=${biomarker/ /_}
    output_file="${folder}predict_${biomarker_str}.pdf"
    prediction/evaluate_predictions.py --predict_biomarker "${biomarker}" --plot_file "${output_file}"
  done
fi

################################ 
###       IMPI               ###
################################
folder="/vol/medic01/users/aschmidt/Dropbox/Documents/Einreichungen/IPMI 2015/figs/"
#folder="/Users/aschmiri/Dropbox/Documents/Einreichungen/IPMI 2015/figs/"

#
# Model combination
#
if [ $experiment -eq "9" ] ; then
  biomarker="MMSE"
  training/plot_models.py -b "${biomarker}" -p cnmci --no_densities --no_sample_lines --ylim 5 35 --plot_file "${folder}model_mmse_cnmci.pdf"
  training/plot_models.py -b "${biomarker}" -p mciad --no_densities --no_sample_lines --ylim 5 35 --plot_file "${folder}model_mmse_mciad.pdf"
  training/plot_models.py -b "${biomarker}" -p joint --no_densities --no_sample_lines --ylim 5 35 --plot_file "${folder}model_mmse_joint.pdf"
fi

#
# Model examples
#
if [ $experiment -eq "10" ] ; then
  for biomarker in "Right Hippocampus" "CDRSB" "MMSE" "D1" "D2"
  do
    biomarker_str=${biomarker/ /_}
    output_file="${folder}model_${biomarker_str}.pdf"
    training/plot_models.py -b "${biomarker}" -p joint --plot_file "${output_file}"
  done
fi

#
# Progression estimates
#
if [ $experiment -eq "11" ] ; then
  for method in  hcv #cog ml vol img all
  do
    i=1
    for visits in "bl" "bl m12" "bl m12 m24"
    do
      output_file="${folder}eval_fitting_${method}_${i}.pdf"
      estimation/estimate_progressions.py ${visits} -m ${method} -p joint --consistent_data --no_analysis --plot_file "${output_file}"
      i=`expr $i + 1`   
    done
  done
fi

