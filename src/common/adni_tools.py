#! /usr/bin/env python
# print __doc__

import os.path
import numpy as np

################################################################################
#
# global paths and variables
#
################################################################################
adni_folder    = '/vol/medic01/users/aschmidt/projects/Data/ADNI'
data_folder    = os.path.join( adni_folder, 'data' )
merge_folder   = os.path.join( adni_folder, 'ADNIMERGE' )
project_folder = '/vol/medic01/users/aschmidt/projects/AgeingAtlas'
param_folder   = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/parameters'
mni_folder     = '/vol/medic01/users/aschmidt/projects/Data/mni_icbm152_nlin_asym_09a'
mni_atlas      = os.path.join( mni_folder, 'mni_icbm152_t1_tal_nlin_asym_09a_brain_scaled.nii' )

################################################################################
#
# make_dir
#
################################################################################
def make_dir( dir1, dir2 ):
    directory = os.path.join( dir1, dir2 )
    if not os.path.exists( directory ):
        os.makedirs( directory )
    return directory

################################################################################
#
# read_list
#
################################################################################
def read_list( query_list, folder ):
    import csv
    files = []
    rids = []
    with open( query_list, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            filename = os.path.join( folder, row[headers.index('Files')] )
            rid = int( row[headers.index('RID')] )
            if os.path.isfile( filename ) and qc_check( rid ):
                files.append( filename )
                rids.append( rid )
    return files, rids

################################################################################
#
# read_list_for_diagnosis
#
################################################################################
def read_list_for_diagnosis( query_list, folder, diagnosis ):
    import csv
    files = []
    rids = []
    with open( query_list, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            filename = os.path.join( folder, row[headers.index('Files')] )
            filename = filename.replace( '.nii.gz', '.dof.gz' )
            rid = int( row[headers.index('RID')] )
            dx = row[headers.index('DX.bl')]
            if diagnosis == 'ALL' or diagnosis in dx or dx in diagnosis:
                if os.path.isfile( filename ):
                    files.append( filename )
                    rids.append( rid )
    return files, rids

################################################################################
#
# read_list_all_data
#
################################################################################
def read_list_all_data( query_list, folder, diagnosis='ALL' ):
    import csv
    files = []
    rids = []
    diagnoses = []
    ages = []
    mmses = []
    
    with open( query_list, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            filename = os.path.join( folder, row[headers.index('Files')] )
            if not os.path.isfile( filename ):
                filename = filename.replace( '.nii.gz', '.dof.gz' )
            if os.path.isfile( filename ):
                rid = row[headers.index('RID')]
                dx = row[headers.index('DX.bl')]
                age = row[headers.index('AGE')]
                mmse = row[headers.index('MMSE')]
                
                if diagnosis == 'ALL' or diagnosis in dx or dx in diagnosis:
                    files.append( filename )
                    rids.append( rid )
                    diagnoses.append( dx )
                    ages.append( age )
                    mmses.append( mmse )
                
    return files, rids, diagnoses, ages, mmses

################################################################################
#
# qc_check
#
################################################################################
def qc_check( rid ):
    excluded_rids = []
    #excluded_rids = [1268,1286,1380,187,359,384,478,513,544,689,741,973]
    if rid in excluded_rids:
        return False
    else:
        return True

################################################################################
#
# check_mask
#
################################################################################
def check_mask( mask ):
    if os.path.isfile( mask ):
        return mask 
    else:
        # Check if ADNIGO mask exists
        mask_GO = mask.replace( 'ADNI2', 'ADNIGO' )
        if os.path.isfile( mask_GO ):
            return mask_GO
        
        # Check if other image with same image ID exists
        import glob
        mask_folder = os.path.dirname( mask )
        masks = glob.glob( mask_folder + '/*' + mask[mask.rfind( '_' ):] )
        if len( masks ) == 0:
            print 'No alternative mask found!'
            return 'none'
        if len( masks ) > 1:
            print 'Multiple alternative mask found! Selecting fist one...'
        mask = masks[0]
        print 'Alternative mask found: ' + mask
        return mask

################################################################################
#
# get_baseline
#
################################################################################
def get_baseline( baseline_folder, study, field_strength, viscode = 'bl' ):

    baseline_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_' + study + '_' + viscode + '_' + field_strength + 'T.csv'
    baseline_files, _ = read_list( baseline_list, baseline_folder )
                
    return np.array( baseline_files )

################################################################################
#
# get_baseline_and_followup
#
################################################################################
def get_baseline_and_followup( baseline_folder, followup_folder, study, field_strength, viscode ):
    baseline_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_' + study + '_bl_' + field_strength + 'T.csv'
    followup_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_' + study + '_' + viscode + '_' + field_strength + 'T.csv'

    baseline_files_unsorted, baseline_rids_unsorted = read_list( baseline_list, baseline_folder )
    followup_files_unsorted, followup_rids_unsorted = read_list( followup_list, followup_folder )

    # For ADNI2 follow-ups, corresponding baselines may be in the ADNIGO folder
    if study == 'ADNI2':
        baseline_list_GO = baseline_list.replace( 'ADNI2', 'ADNIGO' )
        baseline_folder_GO = baseline_folder.replace( 'ADNI2', 'ADNIGO' )
        
        baseline_files_GO, baseline_rids_GO = read_list( baseline_list_GO, baseline_folder_GO )
        
        baseline_files_unsorted += baseline_files_GO
        baseline_rids_unsorted += baseline_rids_GO
                        
    baseline_files = []
    followup_files = []
    
    for i in range( len(baseline_files_unsorted) ):
        baseline_rid = baseline_rids_unsorted[i]
        for j in range( len(followup_files_unsorted) ):
            followup_rid = followup_rids_unsorted[j]
            if baseline_rid == followup_rid:
                baseline_files.append( baseline_files_unsorted[i] )
                followup_files.append( followup_files_unsorted[j] )
                
    return np.array( baseline_files), np.array( followup_files )

################################################################################
#
# get_baseline_transformations
#
################################################################################
def get_baseline_transformations( dof_folder_adni1, dof_folder_adni2, viscode, diagnosis='ALL' ):
    adni1_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI1_' + viscode + '_1.5T.csv'
    adni2_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI2_' + viscode + '_3T.csv'

    velocities_ADNI1, _ = read_list_for_diagnosis( adni1_list, dof_folder_adni1, diagnosis )
    velocities_ADNI2, _ = read_list_for_diagnosis( adni2_list, dof_folder_adni2, diagnosis )
    
    print 'Found ' + str(len( velocities_ADNI1 )) + ' velocities in ADNI1...'
    print 'Found ' + str(len( velocities_ADNI2 )) + ' velocities in ADNI2/GO...'
    
    return np.array( velocities_ADNI1 + velocities_ADNI2 )

################################################################################
#
# get_baseline_transformations_and_rids
#
################################################################################
def get_baseline_transformations_and_rids( dof_folder_adni1, dof_folder_adni2, viscode, diagnosis='ALL' ):
    adni1_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI1_' + viscode + '_1.5T.csv'
    adni2_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI2_' + viscode + '_3T.csv'

    velocities_ADNI1, rids_ADNI1 = read_list_for_diagnosis( adni1_list, dof_folder_adni1, diagnosis )
    velocities_ADNI2, rids_ADNI2 = read_list_for_diagnosis( adni2_list, dof_folder_adni2, diagnosis )
    
    print 'Found ' + str(len( velocities_ADNI1 )) + ' velocities in ADNI1...'
    print 'Found ' + str(len( velocities_ADNI2 )) + ' velocities in ADNI2/GO...'
    
    velocities = np.array( velocities_ADNI1 + velocities_ADNI2 )
    rids = np.array( rids_ADNI1 + rids_ADNI2, dtype='int' )
    
    return velocities, rids

################################################################################
#
# get_all_data
#
################################################################################
def get_all_data( data_folder_adni1, data_folder_adni2, data_folder_adniG, viscode, diagnosis='ALL' ):
    adni1_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI1_' + viscode + '_1.5T.csv'
    adni2_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI2_' + viscode + '_3T.csv'
    adniG_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNIGO_' + viscode + '_3T.csv'

    files1, rid1, diagnoses1, age1, mmse1 = read_list_all_data( adni1_list, data_folder_adni1, diagnosis )
    files2, rid2, diagnoses2, age2, mmse2 = read_list_all_data( adni2_list, data_folder_adni2, diagnosis )
    if viscode == 'bl':
        filesG, ridG, diagnosesG, ageG, mmseG = read_list_all_data( adniG_list, data_folder_adniG, diagnosis )
    else:
        filesG = []; ridG = []; diagnosesG = []; ageG = []; mmseG = []
        
    print 'Found ' + str(len( files1 )) + ' files in ADNI1...'
    print 'Found ' + str(len( files2 )) + ' files in ADNI2...'
    print 'Found ' + str(len( filesG )) + ' files in ADNIGO...'
        
    files = np.array( files1 + files2 + filesG )
    rids = np.array( rid1 + rid2 + ridG, dtype='int' )
    diagnoses = np.array( diagnoses1 + diagnoses2 + diagnosesG )
    ages = np.array( age1 + age2 + ageG, dtype='float' )
    mmses = np.array( mmse1 + mmse2 + mmseG, dtype='int' )
    
    return files, rids, diagnoses, ages, mmses 

################################################################################
#
# find_images_with_dof
#
################################################################################
def find_images_with_dof( image_files, dof_folder ):
    image_files_with_dof = []
    dof_files = []
    for i in range( len(image_files) ):
        image_base = os.path.basename( image_files[i] )
        dof = os.path.join( dof_folder, image_base )
        dof = dof.replace('.nii.gz', '.dof.gz')
        if os.path.isfile( dof ):
            image_files_with_dof.append( image_files[i] )
            dof_files.append( dof )
            
    return image_files_with_dof, dof_files


################################################################################
#
# get_years_after_symptoms
#
################################################################################
def get_years_after_symptoms():
    import csv
    import sqlite3
    
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    
    #
    # Read dates database 
    #
    cur.execute("CREATE TABLE dates (rid INT PRIMARY KEY, \
                                     viscode TEXT, \
                                     year_symptom INT, \
                                     year_mci INT, \
                                     year_ad INT)")

    file_dates = os.path.join( project_folder, 'lists/ADNI_diag_dates.csv' )
    
    with open( file_dates, 'rb' ) as csvfile:
        entries = csv.DictReader( csvfile )
        for entry in entries:
            rid          = entry['Subject ID']
            viscode      = entry['Visit ID']
            year_symptom = entry['Year cognitive symptoms began']
            year_mci     = entry['Year MCI diagnosed']
            year_ad      = entry['Year AD diagnosed']
            if viscode == 'sc' or viscode == 'v01':
                cur.execute("INSERT INTO dates VALUES (?,?,?,?,?)", (rid, viscode, year_symptom, year_mci, year_ad) )
                
    cur.execute("UPDATE dates SET year_symptom=NULL WHERE year_symptom=''" )
    cur.execute("UPDATE dates SET year_mci=NULL WHERE year_mci=''" )
    cur.execute("UPDATE dates SET year_ad=NULL WHERE year_ad=''" )
        
    #
    # Read dpi database 
    #
    cur.execute("CREATE TABLE dpis (rid INT PRIMARY KEY, \
                                    dpi FLOAT)")

    file_dates = os.path.join( project_folder, 'atlas/model_1/data_sym_m24_AD.csv' )
    
    with open( file_dates, 'rb' ) as csvfile:
        entries = csv.DictReader( csvfile )
        for entry in entries:
            rid      = entry['RID']
            dpi      = entry['DPI']
            cur.execute("INSERT INTO dpis VALUES (?,?)", (rid, dpi) )
    
    #
    # Read summary database 
    #
    cur.execute("CREATE TABLE summary (rid INT PRIMARY KEY, \
                                       filename TEXT, \
                                       diagnosis TEXT, \
                                       scandate TEXT, \
                                       age FLOAT, \
                                       apoe4 INT, \
                                       cdrsb FLOAT, \
                                       adas11 FLOAT, \
                                       adas13 FLOAT, \
                                       faq INT, \
                                       mmse INT, \
                                       moca INT)")

    def add_summary( folder, filename ):
        with open( filename, 'rb' ) as csvfile:
            entries = csv.DictReader( csvfile )
            for entry in entries:
                filename = os.path.join( folder, entry['Files'] )
                if not os.path.isfile( filename ):
                    filename = filename.replace( '.nii.gz', '.dof.gz' )
                if os.path.isfile( filename ):
                    rid      = entry['RID']
                    dx       = entry['DX.bl']
                    scandate = entry['ScanDate']
                    age      = entry['AGE']
                    apoe4    = entry['APOE4']
                    cdrsb    = entry['CDRSB']
                    adas11   = entry['ADAS11']
                    adas13   = entry['ADAS13']
                    faq      = entry['FAQ']
                    mmse     = entry['MMSE']
                    moca     = entry['MOCA']
                    cur.execute("INSERT INTO summary VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                                 (rid, filename, dx, scandate, age, apoe4, cdrsb, adas11, adas13, faq, mmse, moca) )

    add_summary( '/vol/medic01/users/aschmidt/projects/Data/ADNI/data/ADNI1/native/images', os.path.join( project_folder, 'lists/query_ADNI1.csv' ) )
    add_summary( '/vol/medic01/users/aschmidt/projects/Data/ADNI/data/ADNI2/native/images', os.path.join( project_folder, 'lists/query_ADNI2.csv' ) )
    add_summary( '/vol/medic01/users/aschmidt/projects/Data/ADNI/data/ADNIGO/native/images', os.path.join( project_folder, 'lists/query_ADNIGO.csv' ) )
    
    cur.execute("UPDATE summary SET apoe4=NULL WHERE apoe4='NA'" )
    cur.execute("UPDATE summary SET cdrsb=NULL WHERE cdrsb='NA'" )
    cur.execute("UPDATE summary SET adas11=NULL WHERE adas11='NA'" )
    cur.execute("UPDATE summary SET adas13=NULL WHERE adas13='NA'" )
    cur.execute("UPDATE summary SET faq=NULL WHERE faq='NA'" )
    cur.execute("UPDATE summary SET mmse=NULL WHERE mmse='NA'" )
    cur.execute("UPDATE summary SET moca=NULL WHERE moca='NA'" )
    
    cur.execute("SELECT rid, scandate, year_symptom, dpi, cdrsb, adas11, adas13, faq, mmse, moca FROM summary NATURAL JOIN dates NATURAL JOIN dpis")
    
    rows = cur.fetchall()
    years  = []
    dpi    = []
    cdrsb  = []
    adas11 = []
    adas13 = []
    faq    = []
    mmse   = []
    moca   = [] 
    for row in rows:
        if row[2] != None:
            print row
            scan_year  = int(row[1][0:4])
            years.append( scan_year - row[2] )
            dpi.append( row[3] )
            cdrsb.append( row[4] )
            adas11.append( row[5] )
            adas13.append( row[6] )
            faq.append( row[7] )
            mmse.append( row[8] )
            moca.append( row[9] )
    return years, dpi, cdrsb, adas11, adas13, faq, mmse, moca
            
        


################################################################################
#
# print_data_for_rid
#
################################################################################   
def print_data_for_rid( rid = 43 ): 
    import rpy2.robjects as robjects
    
    # Lead database
    robjects.r['load']( os.path.join( merge_folder, 'data/adnimerge.rdata' ) )
    adnimerge = robjects.r['adnimerge']
    
    # Print header
    for col in range(adnimerge.ncol):
        print adnimerge.colnames[col],
    print ''
    
    # Print each row
    rids = adnimerge[adnimerge.colnames.index('RID')]
    for row in range(adnimerge.nrow):
        if rids[row] == rid:
            for col in range(adnimerge.ncol):
                print adnimerge[col][row],
            print ''


################################################################################
#
# main
#
################################################################################
if __name__ == "__main__":
    get_years_after_symptoms()
    #print_data_for_rid()
#     baseline_folder = '/vol/vipdata/data/ADNI/data/ADNI1/native/images_unstripped/'
#     followup_folder = '/vol/vipdata/data/ADNI/data/ADNI1/native/images_unstripped/'
#  
#     [baseline_files, followup_files] = get_baseline_and_followup( baseline_folder, followup_folder, 'ADNI1', '1.5' )
#     print 'Found ' + str(len(baseline_files)) + ' image pairs in ADNI1 (1.5T)'
#  
#     #[baseline_files, followup_files] = get_baseline_and_followup( baseline_folder, followup_folder, 'ADNI1', '3' )
#     #print 'Found ' + str(len(baseline_files)) + ' image pairs in ADNI1 (3T)'
#      
#     baseline_folder = '/vol/vipdata/data/ADNI/data/ADNI2/native/images_unstripped/'
#     followup_folder = '/vol/vipdata/data/ADNI/data/ADNI2/native/images_unstripped/'
#  
#     [baseline_files, followup_files] = get_baseline_and_followup( baseline_folder, followup_folder, 'ADNI2', '3' )
#     print 'Found ' + str(len(baseline_files)) + ' image pairs in ADNI2'
