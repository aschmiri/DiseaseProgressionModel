#! /usr/bin/env python
# print __doc__

import sys
import os.path
import csv
    
if len( sys.argv ) < 3:
    print 'Usage: applyBrainMasks.py <study> <field_strength>'
    exit()
  
study = sys.argv[1]
fs = sys.argv[2]

listt = '/vol/biomedic/users/aschmidt/AgingAtlas/lists/query_' + study + '_m24_' + fs + 'T.csv'
folder = '/vol/biomedic/users/aschmidt/ADNI/data/' + study + '/baseline_linear_niftyreg/images_unstripped'

# Read CSV file
num_files = 0
found = []
with open( listt, 'rb' ) as csvfile:
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
