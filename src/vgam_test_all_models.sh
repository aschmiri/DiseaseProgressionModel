#! /bin/bash
#method=$1

#for method in cog mbl long img all
#do
#  for visits in  "bl" "m12" "m24" "bl m12" "m12 m24" "bl m12 m24"
#  do
#    fitting/vgam_estimate_progress.py ${visits} -m ${method} --no_plot
#  done
#done

for biomarker in "MMSE" "FAQ" "CDRSB" "Right Hippocampus"
do
  for method in cog mbl long img all
  do
    for visits in  "bl" "m12" "m24" "bl m12" "m12 m24" "bl m12 m24"
    do
      fitting/vgam_predict_biomarker_value.py ${visits} -m ${method} -p "${biomarker}" --consistent_data --exclude_cn --latex_file results.tex --no_plot
    done
  done
done

#for method in cog mbl all img long
#do
#  for visits in "bl m12" "bl m12 m24" "m12 m24"
#  do
#    fitting/vgam_estimate_progress.py ${visits} -m ${method} --estimate_dpr --no_plot
#  done
#done

for biomarker in "MMSE" "FAQ" "CDRSB" "Right Hippocampus"
do
  for method in cog mbl all img long
  do
    for visits in "bl m12" "bl m12 m24" "m12 m24"
    do
      fitting/vgam_predict_biomarker_value.py ${visits} -m ${method} -p "${biomarker}" --estimate_dpr --consistent_data --exclude_cn --latex_file results.tex --no_plot
    done
  done
done

