#! /usr/bin/env python
# print __doc__
import csv
import numpy as np
import matplotlib.pyplot as plt


def main():
    datafile_m0 = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_0/data_sym_m24_AD.csv'
    datafile_m1 = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_1/data_sym_m24_AD.csv'
    datafile_m2 = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_2/data_sym_m24_AD.csv'

    #
    # Read data
    states_m0_cn = read_data(datafile_m0, ['CN'])
    states_m0_mci = read_data(datafile_m0, ['EMCI', 'LMCI'])
    states_m0_ad = read_data(datafile_m0, ['AD'])

    states_m1_cn = read_data(datafile_m1, ['CN'])
    states_m1_mci = read_data(datafile_m1, ['EMCI', 'LMCI'])
    states_m1_ad = read_data(datafile_m1, ['AD'])

    # states_m2_cn  = read_data( datafile_m2, ['CN'] )
    # states_m2_mci  = read_data( datafile_m2, ['EMCI', 'LMCI'] )
    # states_m2_ad  = read_data( datafile_m2, ['AD'] )

    plt.scatter(states_m0_ad, states_m1_ad, color=(0, 0.38, 0.48), alpha=0.2)

    coefs_m1 = np.polyfit(states_m0_ad, states_m1_ad, 2)
    # coefs_m2 =  np.polyfit( states_m0_ad, states_m2_ad, 2)
    print coefs_m1  # , '  ', coefs_m2

    x = np.arange(0, 16)
    y = np.square(x) * coefs_m1[0] + x * coefs_m1[1] + coefs_m1[2]
    plt.plot(x, y, color=(0, 0.38, 0.48))

    # y = np.square(x) * coefs_m2[0] + x * coefs_m2[1] + coefs_m2[2]
    # plt.plot( x, y, '--', color=(0,0.38,0.48) )

    plt.plot([-0.5, 15.5], [-0.5, 15.5], '--', color='grey')

    plt.xlabel('Virtual disease state, iteration 0')
    plt.ylabel('Virtual disease state, iteration 1')
    plt.xlim([-0.5, 15.5])
    plt.ylim([-0.5, 15.5])

    plt.show()


def read_data(datafile, diagnoses):
    states = []

    with open(datafile, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        headers = reader.next()
        for row in reader:
            dx = row[headers.index('DX.bl')]
            if dx in diagnoses:
                states.append(float(row[headers.index('DPI')]))

    return states


if __name__ == '__main__':
    main()
