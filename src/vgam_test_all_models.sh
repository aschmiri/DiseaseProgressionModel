#! /bin/bash
if [ -z $1 ] ; then
  methods=(cog ml vol img all)
else
  methods=($1)
fi

for method in ${methods[*]}
do
  for visits in  "bl" "bl m12" "bl m12 m24" "m12" "m24" "m12 m24"
  do
    estimation/estimate_progressions.py ${visits} -m ${method} --no_plot
  done
done

for biomarker in "MMSE" "FAQ" "CDRSB" "Right Hippocampus"
do
  for method in ${methods[*]}
  do
    for visits in  "bl" "bl m12" "bl m12 m24" "m12" "m24" "m12 m24"
    do
      prediction/predict_biomarker_values.py ${visits} -m ${method} -p "${biomarker}" --consistent_data --exclude_cn --no_plot
    done
  done
done

for method in ${methods[*]}
do
  for visits in "bl m12 m24" "bl m12" "m12 m24"
  do
    estimation/estimate_progressions.py ${visits} -m ${method} --estimate_dpr --no_plot
  done
done

for biomarker in "MMSE" "FAQ" "CDRSB" "Right Hippocampus"
do
  for method in ${methods[*]}
  do
    for visits in "bl m12 m24" "bl m12" "m12 m24"
    do
      prediction/predict_biomarker_values.py ${visits} -m ${method} -p "${biomarker}" --estimate_dpr --consistent_data --exclude_cn --no_plot
    done
  done
done

