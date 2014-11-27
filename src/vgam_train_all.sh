#! /bin/bash
if [ -z $1 ] ; then
  phases=(cnmci mciad joint)
else
  phases=($1)
fi

for phase in ${phases[*]}
do
  training/train_models.py -p ${phase}
  training/evaluate_models.py -p ${phase}
  training/plot_models.py -p ${phase} --save_plot
done

