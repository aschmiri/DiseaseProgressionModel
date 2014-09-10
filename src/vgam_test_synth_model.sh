#! /bin/bash

for biomarker in synth_e1 synth_e2
do
  for num_subjects in 1000 # 250 500 1000 1500 2000
  do
    fitting/vgam_generate_synth_data.py -b ${biomarker} -n ${num_subjects} #--uniform_progression
    fitting/vgam_estimate_curves.py synth -b ${biomarker}
    fitting/vgam_plot_curves.py synth -b ${biomarker} --plot_synth_model --no_extrapolation
  done
done
