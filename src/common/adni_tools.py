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

query_list     = os.path.join( project_folder, 'lists/query_ADNI.csv' )

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
def read_list( folder, diagnosis = 'ALL', study = 'ALL', viscode = 'ALL' ):
    files, rids, _, _, _ = read_list_all_data( folder, diagnosis, study, viscode )
    return files, rids

################################################################################
#
# read_list_all_data
#
################################################################################
def read_list_all_data( folder, diagnosis='ALL', study='ALL', viscode='ALL' ):
    import csv
    import re
    files = []
    rids = []
    diagnoses = []
    ages = []
    mmses = []
    
    with open( query_list, 'rb' ) as csvfile:
        reader = csv.reader( csvfile, delimiter=',' )
        headers = reader.next()
        for row in reader:
            vc = row[headers.index('VISCODE')]
            fs = row[headers.index('MagStrength')]
            std = row[headers.index('COLPROT')]
            rid = row[headers.index('RID')]
            dx = row[headers.index('DX.bl')]
            age = row[headers.index('AGE.scan')]
            mmse = row[headers.index('MMSE')]
                
            if ((diagnosis == 'ALL' or diagnosis in dx or dx in diagnosis) and
                (study == 'ALL' or study == std) and
                ((std == 'ADNI1' and float(fs) < 2.0) or (std != 'ADNI1' and float(fs) > 2.0)) and
                (viscode == 'ALL' or viscode == vc or (viscode == 'fu' and re.match('m[0-9][0-9]', vc)))):
                filename_base = os.path.basename( row[headers.index('Files')] )
                filename = os.path.join( folder, filename_base )
                if not os.path.isfile( filename ):
                    filename = filename.replace( '.nii.gz', '.dof.gz' )
                if os.path.isfile( filename ):
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
            return None
        
        if len( masks ) > 1:
            print 'WARNING: Multiple alternative mask found! Selecting fist one...'

        mask = masks[0]
        print 'Alternative mask found: ' + mask
        return mask

################################################################################
#
# get_images
#
################################################################################
def get_images( folder, study, viscode ):
    baseline_files, _ = read_list( folder, study=study, viscode=viscode )
                
    return np.array( baseline_files )

################################################################################
#
# get_baseline
#
################################################################################
def get_baseline( baseline_folder, study ):
    baseline_files = get_images( baseline_folder, study=study, viscode='bl' )
                
    return np.array( baseline_files )

################################################################################
#
# get_followup
#
################################################################################
def get_followup( baseline_folder, study, viscode = 'fu' ):
    followup_files = get_images( baseline_folder, study=study, viscode=viscode )
                
    return np.array( followup_files )

################################################################################
#
# get_baseline_and_followup
#
################################################################################
def get_baseline_and_followup( baseline_folder, followup_folder, study, viscode ):
    baseline_files_unsorted, baseline_rids_unsorted = read_list( baseline_folder, study=study, viscode='bl' )
    followup_files_unsorted, followup_rids_unsorted = read_list( followup_folder, study=study, viscode=viscode )
    
    # For ADNI2 follow-ups, corresponding baselines may be in the ADNIGO folder
    if study == 'ADNI2':
        baseline_folder_GO = baseline_folder.replace( 'ADNI2', 'ADNIGO' )
        baseline_files_GO, baseline_rids_GO = read_list( baseline_folder_GO, study='ADNIGO', viscode='bl' )
        
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
def get_baseline_transformations( dof_folder_adni1, dof_folder_adni2, viscode, diagnosis = 'ALL' ):
    velocities_ADNI1, _ = read_list( dof_folder_adni1, diagnosis=diagnosis, study='ADNI1', viscode=viscode )
    velocities_ADNI2, _ = read_list( dof_folder_adni2, diagnosis=diagnosis, study='ADNI2', viscode=viscode )
    
    print 'Found', len( velocities_ADNI1 ), 'velocities in ADNI1...'
    print 'Found', len( velocities_ADNI2 ), 'velocities in ADNI2/GO...'
    
    return np.array( velocities_ADNI1 + velocities_ADNI2 )

################################################################################
#
# get_baseline_transformations_and_rids
#
################################################################################
def get_baseline_transformations_and_rids( dof_folder_adni1, dof_folder_adni2, viscode, diagnosis='ALL' ):
    velocities_ADNI1, rids_ADNI1 = read_list( dof_folder_adni1, diagnosis=diagnosis, study='ADNI1', viscode=viscode )
    velocities_ADNI2, rids_ADNI2 = read_list( dof_folder_adni2, diagnosis=diagnosis, study='ADNI2', viscode=viscode )
    
    print 'Found', len( velocities_ADNI1 ), 'velocities in ADNI1...'
    print 'Found', len( velocities_ADNI2 ), 'velocities in ADNI2/GO...'
    
    velocities = np.array( velocities_ADNI1 + velocities_ADNI2 )
    rids = np.array( rids_ADNI1 + rids_ADNI2, dtype='int' )
    
    return velocities, rids

################################################################################
#
# get_all_data
#
################################################################################
def get_all_data( data_folder_adni1, data_folder_adni2, data_folder_adniG, viscode, diagnosis='ALL' ):
    files1, rid1, diagnoses1, age1, mmse1 = read_list_all_data( data_folder_adni1, diagnosis=diagnosis, study='ADNI1', viscode=viscode  )
    files2, rid2, diagnoses2, age2, mmse2 = read_list_all_data( data_folder_adni2, diagnosis=diagnosis, study='ADNI2', viscode=viscode  )
    if viscode == 'bl':
        filesG, ridG, diagnosesG, ageG, mmseG = read_list_all_data( data_folder_adniG, diagnosis=diagnosis, study='ADNIGO', viscode=viscode  )
    else:
        filesG = []; ridG = []; diagnosesG = []; ageG = []; mmseG = []
        
    print 'Found', len( files1 ), 'files in ADNI1...'
    print 'Found', len( files2 ), 'files in ADNI2...'
    print 'Found', len( filesG ), 'files in ADNIGO...'
        
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
# check_image_data
#
################################################################################ 
def check_image_data():
    studies = ['ALL', 'ADNI1', 'ADNI2', 'ADNIGO']
    followup_viscodes = ['m06', 'm12', 'm18', 'm24', 'm36', 'm48', 'm60'] 
    
    for study in studies:
        folder = os.path.join( data_folder, study, 'native/images_unstripped' )
        
        image_files = get_images( folder, study=study, viscode='ALL' )
        print 'Found', len(image_files), 'images in', study
        
        baseline_files = get_baseline( folder, study=study )
        print 'Found', len(baseline_files), 'baseline images in', study
    
        followup_files = get_followup( folder, study=study )
        print 'Found', len(followup_files), 'followup images in', study
        
        for vc in followup_viscodes:
            followup_files = get_followup( folder, study = study, viscode=vc )
            print 'Found', len(followup_files), 'followup images in', study, '(' + vc + ')'
            
        baseline_files, _ = get_baseline_and_followup( folder, folder, study=study, viscode='fu' )
        print 'Found', len(baseline_files), 'image pairs in', study, '(fu)'


################################################################################
#
# main
#
################################################################################
if __name__ == "__main__":
    check_image_data()

