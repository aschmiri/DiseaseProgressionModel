#! /bin/bash
method=$1

for method in cog mbl vol img all
do
  for visits in  "bl" "m12" "m24" "bl m12" "m12 m24" "bl m12 m24"
  do
    estimation/estimate_progressions.py ${visits} -m ${method} --no_plot
  done
done

for biomarker in "MMSE" "FAQ" "CDRSB" "Right Hippocampus"
do
  for method in cog mbl vol img all
  do
    for visits in  "bl" "m12" "m24" "bl m12" "m12 m24" "bl m12 m24"
    do
      prediction/predict_biomarker_values.py ${visits} -m ${method} -p "${biomarker}" --consistent_data --exclude_cn --latex_file results.tex --no_plot
    done
  done
done

for method in cog mbl vol img all
do
  for visits in "bl m12" "bl m12 m24" "m12 m24"
  do
    estimation/estimate_progressions.py ${visits} -m ${method} --estimate_dpr --no_plot
  done
done

for biomarker in "MMSE" "FAQ" "CDRSB" "Right Hippocampus"
do
  for method in cog mbl vol img all
  do
    for visits in "bl m12" "bl m12 m24" "m12 m24"
    do
      prediction/predict_biomarker_values.py ${visits} -m ${method} -p "${biomarker}" --estimate_dpr --consistent_data --exclude_cn --latex_file results.tex --no_plot
    done
  done
done

