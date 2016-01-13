#! /bin/bash
experiment=${1-"0"}

################################ 
###           PLOS           ###
################################

# Classify diagnoses
if [ ${experiment} -eq "1" ] ; then
  texfile="./results_diagnosis.tex"
  for classifier in lda # svm lsvm rf
  do
    for method in img vol ml #img
    do
      classification/classify_diagnoses.py bl -m ${method} -p mciad -c ${classifier} --consistent_data --latex_file ${texfile}
      classification/classify_diagnoses.py bl m12 -m ${method} -p mciad -c ${classifier} --consistent_data --latex_file ${texfile}
      classification/classify_diagnoses.py bl m12 m24 -m ${method} -p mciad -c ${classifier} --consistent_data --latex_file ${texfile}
    done
  done
fi


################################ 
###           IPMI           ###
################################

# Classify diagnoses
if [ ${experiment} -eq "10" ] ; then
  texfile="results_diagnosis.tex"
  for classifier in svm lsvm lda rf
  do
    for method in cog ml vol img img2 #all
    do
      classification/classify_diagnoses.py bl -m ${method} -p joint -c ${classifier} --consistent_data --latex_file ${texfile}
      classification/classify_diagnoses.py bl m12 m24 -m ${method} -p joint -c ${classifier} --consistent_data --latex_file ${texfile}
      classification/classify_diagnoses.py bl m12 m24 -m ${method} -p joint -c ${classifier} --estimate_dpr --consistent_data --latex_file ${texfile}
    done
  done
fi

# Classify converters
if [ ${experiment} -eq "11" ] ; then
  texfile="results_convert.tex"
  for method in cog ml vol img img2 #all
  do
    classification/classify_converters.py bl -m ${method} -p joint --consistent_data --latex_file ${texfile}
    classification/classify_converters.py m12 -m ${method} -p joint --consistent_data --latex_file ${texfile}
    classification/classify_converters.py bl m12 -m ${method} -p joint --consistent_data --latex_file ${texfile}
    classification/classify_converters.py bl m12 -m ${method} -p joint --estimate_dpr --consistent_data --latex_file ${texfile}
  done
fi

# Classify rapid decline
if [ ${experiment} -eq "12" ] ; then
  texfile="results_decline.tex"
  for method in all #cog ml vol img img2 all
  do
    classification/classify_rapid_decline.py bl -m ${method} -p joint --consistent_data --latex_file ${texfile}
    classification/classify_rapid_decline.py m12 -m ${method} -p joint --consistent_data --latex_file ${texfile}
    classification/classify_rapid_decline.py bl m12 -m ${method} -p joint --consistent_data --latex_file ${texfile}
    classification/classify_rapid_decline.py bl m12 -m ${method} -p joint --estimate_dpr --consistent_data  --latex_file ${texfile}
  done
fi

