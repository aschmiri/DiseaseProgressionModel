#! /usr/bin/env python2.7
import os.path
import argparse
import csv
import sqlite3
import pickle
from common import adni_tools as adni
from common import log as log


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--estimates_basename', type=str, default='eval_cog_with_rids.p', help='the basename of the file with the DPI estimates')
    parser.add_argument('--output_basename', type=str, default='dpis_for_atlas.csv', help='the basename of the output CSV file')
    args = parser.parse_args()

    # Setup data
    folder_lin = adni.data_folder + '/STUDY/MNI152_linear/images'

    estimates_file = os.path.join(adni.project_folder, 'eval', args.estimates_basename)
    (dpis, _, rids) = pickle.load(open(estimates_file, 'rb'))

    # Setup output file
    output_file = os.path.join(adni.project_folder, 'lists', args.output_basename)
    writer = csv.writer(open(output_file, 'wb'), delimiter=',')
    writer.writerow(['rid', 'dpi', 'filename'])

    # Setup DB
    con = sqlite3.connect(os.path.join(adni.project_folder, 'lists', 'adni.db'))
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    for rid, dpi in zip(rids, dpis):
        cur.execute("SELECT study_bl, filename FROM adnimerge WHERE viscode = 'bl' AND rid = " + str(rid))
        visits = cur.fetchall()
        if len(visits) != 1:
            print log.ERROR, 'Wrong number of visits: ', visits
        else:
            study = visits[0]['study_bl']
            filename_base = visits[0]['filename']
            filename = os.path.join(folder_lin.replace('STUDY', study), filename_base)

            if not os.path.isfile(filename):
                print log.WARNING, 'File not found: {0}'.format(filename)
            else:
                print log.INFO, 'File found: {0}'.format(filename)
                writer.writerow([rid, dpi, filename])


if __name__ == '__main__':
    main()
