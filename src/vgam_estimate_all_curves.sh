#! /bin/bash

process_group() {
    python ./fitting/vgam_estimate_curves.py -m $1 -i $2
    python ./fitting/vgam_evaluate_curves.py -m $1 -i $2
    python ./fitting/vgam_plot_curves.py -m $1 -i $2 --plot_mu --save_file
}


for it in 0 1 2
do
    process_group cog ${it}
    process_group mbl ${it}
    process_group graph ${it}
done

for it in 0 1 2
do
    for bio in long reg cons
    do
        process_group ${bio} ${it}
    done
done
