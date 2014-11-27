#! /bin/bash

# Classify diagnoses
texfile="results_diagnosis.tex"
for classifier in svm lda lsvm svm rf
do
  for method in cog ml vol img img2 #all
  do
    classification/classify_diagnoses.py bl -m ${method} -p joint -c ${classifier} --consistent_data --latex_file ${texfile}
    classification/classify_diagnoses.py bl m12 m24 -m ${method} -p joint -c ${classifier} --consistent_data --latex_file ${texfile}
    classification/classify_diagnoses.py bl m12 m24 -m ${method} -p joint -c ${classifier} --consistent_data --estimate_dpr --latex_file ${texfile}
  done
done

# Classify converters
texfile="results_convert.tex"
for classifier in svm lda lsvm svm rf
do
  for method in cog ml vol img img2 #all
  do
    classification/classify_converters.py bl -m ${method} -p joint -c ${classifier} --consistent_data --latex_file ${texfile}
  done
done

# Classify rapid decline
texfile="results_decline.tex"
for classifier in svm lda lsvm svm rf
do
  for method in cog ml vol img img2 #all
  do
    classification/classify_rapid_decline.py bl -m ${method} -p joint -c ${classifier} --consistent_data --latex_file ${texfile}
  done
done
