#! /usr/bin/env python2.7
import csv
import numpy as np
import matplotlib.pyplot as plt


def main():
    datafile_bl = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/atlas/model_0/data_m24_MCI.csv'

    #
    # Read data
    states_emci = read_data(datafile_bl, ['EMCI'])
    states_lmci = read_data(datafile_bl, ['LMCI'])

    x = np.arange(-0.5, 11, 1)

    plt.hist(states_emci, x, normed=True, color='green', alpha=0.7)
    plt.hist(states_lmci, x, normed=True, color='red', alpha=0.7)

    plt.xlabel('Virtual disease state')
    plt.ylabel('Probability')
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
