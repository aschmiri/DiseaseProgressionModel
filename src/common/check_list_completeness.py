#! /usr/bin/env python
# print __doc__

import argparse
import os.path
import csv
    
parser = argparse.ArgumentParser()
parser.add_argument( 'study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO' )
parser.add_argument( 'field_strength', type=str,  help='the field strength, usually 1.5 for ADNI1 and 3 otherwise' )
a = parser.parse_args()

data_list = '/vol/biomedic/users/aschmidt/AgingAtlas/lists/query_' + a.study + '_m24_' + a.field_strength + 'T.csv'
folder = '/vol/biomedic/users/aschmidt/ADNI/data/' + a.study + '/baseline_linear/images_unstripped'

# Read CSV file
num_files = 0
found = []
with open( data_list, 'rb' ) as csvfile:
    reader = csv.reader( csvfile, delimiter=',' )
    headers = reader.next()
    for row in reader:
        filename = os.path.join( folder, row[headers.index('Files')] )
        if not os.path.isfile( filename ):
            print 'File not found: ' + filename
        else:
            if filename in found:
                print 'Double entry:   ' + filename
            else:
                found.append(filename)
                num_files += 1

print str(num_files)
