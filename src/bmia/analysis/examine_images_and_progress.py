#! /usr/bin/env python2.7
import os.path
import csv
import numpy as np
from subprocess import call
from bmia.common import adni_tools as adni

EXEC_RVIEW = 'rview'


def main():
    datafile_bl = os.path.join(adni.project_folder, 'atlas/data_bl.csv')
    image_folder = os.path.join(adni.data_folder, 'ADNI1/MNI152_linear/images')

    def read_data(datafile, diagnoses):
        images = []
        virtps = []

        with open(datafile, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            headers = reader.next()
            for row in reader:
                dx = row[headers.index('DX.bl')]
                dof = row[headers.index('FILE')]
                basename = os.path.basename(dof).replace('.dof.gz', '.nii.gz')
                image = os.path.join(image_folder, basename)
                if dx in diagnoses and os.path.exists(image):
                    images.append(image)
                    virtps.append(float(row[headers.index('VIRT_P')]))

        return (images, virtps)

    images, virtp = read_data(datafile_bl, ['AD'])
    images = np.array(images)
    virtp = np.array(virtp)

    indices = np.argsort(virtp)
    indices = indices[::-1]

    images_sorted = images[indices]
    virtp_sorted = virtp[indices]

    for i in range(len(images_sorted)):
        print 'VIRT_P: ' + str(virtp_sorted[i]) + ' ' + images_sorted[i]
        call([EXEC_RVIEW, images_sorted[i]])


if __name__ == '__main__':
    main()
