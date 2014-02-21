#! /usr/bin/env python
# print __doc__

import os.path

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
    return (files, rids)

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
            if os.path.isfile( filename ) and (diagnosis == 'ALL' or diagnosis in dx):
                files.append( filename )
                rids.append( rid )
    return files, rids

################################################################################
#
# read_list_all_data
#
################################################################################
def read_list_all_data( query_list, folder ):
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
                diagnosis = row[headers.index('DX.bl')]
                age = row[headers.index('AGE')]
                mmse = row[headers.index('MMSE')]
                
                files.append( filename )
                rids.append( rid )
                diagnoses.append( diagnosis )
                ages.append( age )
                mmses.append( mmse )
                
    return (files, rids, diagnoses, ages, mmses)

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
                
    return baseline_files

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
                
    return (baseline_files, followup_files)

################################################################################
#
# get_velocities
#
################################################################################
def get_velocities( dof_folder_adni1, dof_folder_adni2, viscode, diagnosis='ALL' ):
    adni1_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI1_' + viscode + '_1.5T.csv'
    adni2_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI2_' + viscode + '_3T.csv'

    velocities_ADNI1, _ = read_list_for_diagnosis( adni1_list, dof_folder_adni1, diagnosis )
    velocities_ADNI2, _ = read_list_for_diagnosis( adni2_list, dof_folder_adni2, diagnosis )
    
    print 'Found ' + str(len( velocities_ADNI1 )) + ' velocities in ADNI1...'
    print 'Found ' + str(len( velocities_ADNI2 )) + ' velocities in ADNI2/GO...'
    
    return velocities_ADNI1 + velocities_ADNI2

################################################################################
#
# get_velocities
#
################################################################################
def get_velocities_and_rids( dof_folder_adni1, dof_folder_adni2, viscode, diagnosis ):
    adni1_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI1_' + viscode + '_1.5T.csv'
    adni2_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI2_' + viscode + '_3T.csv'

    velocities_ADNI1, rids_ADNI1 = read_list_for_diagnosis( adni1_list, dof_folder_adni1, diagnosis )
    velocities_ADNI2, rids_ADNI2 = read_list_for_diagnosis( adni2_list, dof_folder_adni2, diagnosis )
    
    print 'Found ' + str(len( velocities_ADNI1 )) + ' velocities in ADNI1...'
    print 'Found ' + str(len( velocities_ADNI2 )) + ' velocities in ADNI2/GO...'
    
    return velocities_ADNI1 + velocities_ADNI2, rids_ADNI1 + rids_ADNI2

################################################################################
#
# get_all_data
#
################################################################################
def get_all_data( data_folder_adni1, data_folder_adni2, data_folder_adniG, viscode ):
    adni1_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI1_' + viscode + '_1.5T.csv'
    adni2_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNI2_' + viscode + '_3T.csv'
    adniG_list = '/vol/medic01/users/aschmidt/projects/AgeingAtlas/lists/query_ADNIGO_' + viscode + '_3T.csv'

    files1, rid1, diagnoses1, age1, mmse1 = read_list_all_data( adni1_list, data_folder_adni1 )
    files2, rid2, diagnoses2, age2, mmse2 = read_list_all_data( adni2_list, data_folder_adni2 )
    if viscode == 'bl':
        filesG, ridG, diagnosesG, ageG, mmseG = read_list_all_data( adniG_list, data_folder_adniG )
    else:
        filesG = []; ridG = []; diagnosesG = []; ageG = []; mmseG = []
        
    print 'Found ' + str(len( files1 )) + ' files in ADNI1...'
    print 'Found ' + str(len( files2 )) + ' files in ADNI2...'
    print 'Found ' + str(len( filesG )) + ' files in ADNIGO...'
    
    return (files1 + files2 + filesG,
            rid1 + rid2 + ridG,
            diagnoses1 + diagnoses2 + diagnosesG,
            age1 + age2 + ageG,
            mmse1 + mmse2 + mmseG)    

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
            
    return (image_files_with_dof, dof_files)

################################################################################
#
# main
#
################################################################################
if __name__ == "__main__":
    baseline_folder = '/vol/vipdata/data/ADNI/data/ADNI1/native/images_unstripped/'
    followup_folder = '/vol/vipdata/data/ADNI/data/ADNI1/native/images_unstripped/'

    [baseline_files, followup_files] = get_baseline_and_followup( baseline_folder, followup_folder, 'ADNI1', '1.5' )
    print 'Found ' + str(len(baseline_files)) + ' image pairs in ADNI1 (1.5T)'

    #[baseline_files, followup_files] = get_baseline_and_followup( baseline_folder, followup_folder, 'ADNI1', '3' )
    #print 'Found ' + str(len(baseline_files)) + ' image pairs in ADNI1 (3T)'
    
    baseline_folder = '/vol/vipdata/data/ADNI/data/ADNI2/native/images_unstripped/'
    followup_folder = '/vol/vipdata/data/ADNI/data/ADNI2/native/images_unstripped/'

    [baseline_files, followup_files] = get_baseline_and_followup( baseline_folder, followup_folder, 'ADNI2', '3' )
    print 'Found ' + str(len(baseline_files)) + ' image pairs in ADNI2'
