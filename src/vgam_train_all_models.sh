#! /bin/bash

process_group() {
    python ./training/train_models.py -m $1
    python ./training/evaluate_models.py -m $1
    python ./training/plot_models.py -m $1 --save_file
}


for bio in cog vol mbl
do
    process_group ${bio} ${it}
done
