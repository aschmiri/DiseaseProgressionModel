#! /bin/bash
if [ -z $1 ] ; then
  methods=(cog ml vol img hcv img2 all)
else
  methods=($1)
fi

for method in ${methods[*]}
do
  for visits in  "bl" "bl m12" "bl m12 m24" #"m12" "m24" "m12 m24"
  do
    estimation/estimate_progressions.py ${visits} -m ${method} -p joint --no_plot --no_analysis
  done
done

for method in ${methods[*]}
do
  for visits in "bl m12" "bl m12 m24" #"m12 m24"
  do
    estimation/estimate_progressions.py ${visits} -m ${method} -p joint --estimate_dpr --no_plot --no_analysis
  done
done

