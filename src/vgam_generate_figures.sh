#! /bin/bash

# Model generation
biomarker="CDRSB"
folder="/vol/medic01/users/aschmidt/Dropbox/Documents/Einreichungen/NeuroImage 2014/figs/"
fitting/vgam_plot_curves.py cog -b "$biomarker" --no_densities --no_extrapolation --no_sample_lines --no_model --output_file="${folder}model_generation_align.pdf"
fitting/vgam_plot_curves.py cog -b "$biomarker" --no_densities --no_extrapolation --no_sample_lines --output_file="${folder}model_generation_fit.pdf"
fitting/vgam_plot_curves.py cog -b "$biomarker" --no_densities --plot_mu --output_file="${folder}model_generation_extrapolate.pdf"
fitting/vgam_plot_curves.py cog -b "$biomarker" --only_densities --output_file="${folder}model_generation_densities.pdf"

# Model examples
for biomarker in "Right Hippocampus" "MMSE" "D1" "D2" "P_D1D2" 
do
  biomarker_str=${biomarker/ /_}
  fitting/vgam_plot_curves.py all -b "$biomarker" --plot_mu --output_file="${folder}model_${biomarker_str}.pdf"
done
