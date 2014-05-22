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

volume_names   = [
    '3rd Ventricle',                  # 0
    '4th Ventricle',                  # 1
    'Right Accumbens Area',           # 2 
    'Left Accumbens Area',            # 3
    'Right Amygdala',                 # 4
    'Left Amygdala',                  # 5
    'Brain Stem',
    'Right Caudate',
    'Left Caudate',
    'Right Cerebellum Exterior',
    'Left Cerebellum Exterior',       #10
    'Right Cerebellum White Matter',
    'Left Cerebellum White Matter',
    'Right Cerebral Exterior',
    'Left Cerebral Exterior',
    'Right Cerebral White Matter',    # 15
    'Left Cerebral White Matter',
    'CSF',
    'Right Hippocampus',
    'Left Hippocampus',
    'Right Inf Lat Vent',             # 20
    'Left Inf Lat Vent',
    'Right Lateral Ventricle',
    'Left Lateral Ventricle',
    'Right Pallidum',
    'Left Pallidum',
    'Right Putamen',
    'Left Putamen',
    'Right Thalamus Proper',
    'Left Thalamus Proper',
    'Right Ventral DC',
    'Left Ventral DC',
    'Right vessel',
    'Left vessel',
    'Optic Chiasm',
    'Cerebellar Vermal Lobules I-V',
    'Cerebellar Vermal Lobules VI-VII',
    'Cerebellar Vermal Lobules VIII-X',
    'Left Basal Forebrain',
    'Right Basal Forebrain',
    'Right ACgG anterior cingulate gyrus',
    'Left ACgG anterior cingulate gyrus',
    'Right AIns anterior insula',
    'Left AIns anterior insula',
    'Right AOrG anterior orbital gyrus',
    'Left AOrG anterior orbital gyrus',
    'Right AnG angular gyrus',
    'Left AnG angular gyrus',
    'Right Calc calcarine cortex',
    'Left Calc calcarine cortex',
    'Right CO central operculum',
    'Left CO central operculum',
    'Right Cun cuneus',
    'Left Cun cuneus',
    'Right Ent entorhinal area',
    'Left Ent entorhinal area',
    'Right FO frontal operculum',
    'Left FO frontal operculum',
    'Right FRP frontal pole',
    'Left FRP frontal pole',
    'Right FuG fusiform gyrus',
    'Left FuG fusiform gyrus',
    'Right GRe gyrus rectus',
    'Left GRe gyrus rectus',
    'Right IOG inferior occipital gyrus',
    'Left IOG inferior occipital gyrus',
    'Right ITG inferior temporal gyrus',
    'Left ITG inferior temporal gyrus',
    'Right LiG lingual gyrus',
    'Left LiG lingual gyrus',
    'Right LOrG lateral orbital gyrus',
    'Left LOrG lateral orbital gyrus',
    'Right MCgG middle cingulate gyrus',
    'Left MCgG middle cingulate gyrus',
    'Right MFC medial frontal cortex',
    'Left MFC medial frontal cortex',
    'Right MFG middle frontal gyrus',
    'Left MFG middle frontal gyrus',
    'Right MOG middle occipital gyrus',
    'Left MOG middle occipital gyrus',
    'Right MOrG medial orbital gyrus',
    'Left MOrG medial orbital gyrus',
    'Right MPoG postcentral gyrus medial segment',
    'Left MPoG postcentral gyrus medial segment',
    'Right MPrG precentral gyrus medial segment',
    'Left MPrG precentral gyrus medial segment',
    'Right MSFG superior frontal gyrus medial segment',
    'Left MSFG superior frontal gyrus medial segment',
    'Right MTG middle temporal gyrus',
    'Left MTG middle temporal gyrus',
    'Right OCP occipital pole',
    'Left OCP occipital pole',
    'Right OFuG occipital fusiform gyrus',
    'Left OFuG occipital fusiform gyrus',
    'Right OpIFG opercular part of the inferior frontal gyrus',
    'Left OpIFG opercular part of the inferior frontal gyrus',
    'Right OrIFG orbital part of the inferior frontal gyrus',
    'Left OrIFG orbital part of the inferior frontal gyrus',
    'Right PCgG posterior cingulate gyrus',
    'Left PCgG posterior cingulate gyrus',
    'Right PCu precuneus',
    'Left PCu precuneus',
    'Right PHG parahippocampal gyrus',
    'Left PHG parahippocampal gyrus',
    'Right PIns posterior insula',
    'Left PIns posterior insula',
    'Right PO parietal operculum',
    'Left PO parietal operculum',
    'Right PoG postcentral gyrus',
    'Left PoG postcentral gyrus',
    'Right POrG posterior orbital gyrus',
    'Left POrG posterior orbital gyrus',
    'Right PP planum polare',
    'Left PP planum polare',
    'Right PrG precentral gyrus',
    'Left PrG precentral gyrus',
    'Right PT planum temporale',
    'Left PT planum temporale',
    'Right SCA subcallosal area',
    'Left SCA subcallosal area',
    'Right SFG superior frontal gyrus',
    'Left SFG superior frontal gyrus',
    'Right SMC supplementary motor cortex',
    'Left SMC supplementary motor cortex',
    'Right SMG supramarginal gyrus',
    'Left SMG supramarginal gyrus',
    'Right SOG superior occipital gyrus',
    'Left SOG superior occipital gyrus',
    'Right SPL superior parietal lobule',
    'Left SPL superior parietal lobule',
    'Right STG superior temporal gyrus',
    'Left STG superior temporal gyrus',
    'Right TMP temporal pole',
    'Left TMP temporal pole',
    'Right TrIFG triangular part of the inferior frontal gyrus',
    'Left TrIFG triangular part of the inferior frontal gyrus',
    'Right TTG transverse temporal gyrus',
    'Left TTG transverse temporal gyrus']

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
    import sqlite3
    con = sqlite3.connect( os.path.join( project_folder, 'lists', 'adni.db' ) )
    #con.row_factory = sqlite3.Row
    cur = con.cursor()
        
    statement_fs = " WHERE ((study IS 'ADNI1' AND fieldstrength < 2.0) OR (study IS NOT 'ADNI1' AND fieldstrength > 2.0))"
    
    if diagnosis == 'ALL':
        statement_diag = ""
    elif diagnosis == 'MCI':
        statement_diag = " AND diagnosis GLOB '*MCI'"
    else:
        statement_diag = " AND diagnosis IS '" + diagnosis + "'"
    
    if study == 'ALL':
        statement_study = ""
    else:
        statement_study = " AND study IS '" + study + "'"
    
    if viscode == 'ALL':
        statement_viscode = ""
    elif viscode == 'fu':
        statement_viscode = " AND viscode GLOB 'm??'"
    else:
        statement_viscode = " AND viscode IS '" + viscode + "'"
    
    cur.execute( "SELECT iid, rid, viscode, study, fieldstrength, age, diagnosis, mmse, filename FROM adnimerge" \
                 + statement_fs + statement_diag + statement_study + statement_viscode )

    files = []
    rids = []
    diagnoses = []
    ages = []
    mmses = []
    
    rows = cur.fetchall()
    for row in rows:
        std = row[3]
        rid = row[1]
        dx = row[6]
        age = row[5]
        mmse = row[7]
               
        if study == 'ALL':
            actual_folder = folder.replace( 'ALL', std )
        else:
            actual_folder = folder
        filename_base = os.path.basename( row[8] )
        filename = os.path.join( actual_folder, filename_base )
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
# find_file
#
################################################################################
def find_file( filename ):
    if os.path.isfile( filename ):
        return filename 
    else:
        # Check if other image with same image ID exists
        import glob
        folder = os.path.dirname( filename )
        files = glob.glob( folder + '/*' + filename[filename.rfind( '_' ):] )
        if len( files ) == 0:
            print 'ERROR: File not found:', filename
            return None
            
        filename = files[0]
        return filename

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
    followup_files_unsorted, followup_rids_unsorted = read_list( followup_folder, study=study, viscode=viscode )
    baseline_files = []
    followup_files = []
    
    import sqlite3
    con = sqlite3.connect( os.path.join( project_folder, 'lists', 'adni.db' ) )
    cur = con.cursor()
    
    for followup_file, followup_rid in zip( followup_files_unsorted, followup_rids_unsorted ):
        cur.execute( "SELECT rid, viscode, study, filename FROM adnimerge WHERE rid = " + str(followup_rid) + " AND viscode = 'bl'" )
        rows = cur.fetchall()
        if len(rows) == 0:
            print 'WARNING: no baseline file found for', followup_file
        elif len(rows) > 1:
            print 'WARNING: multiple baseline files found for', followup_file
        else:
            baseline_file = os.path.join( baseline_folder.replace( study, rows[0][2] ), rows[0][3] )
            if not os.path.isfile( baseline_file ):
                print 'WARNING: baseline file not found', followup_file
            else:
                baseline_files.append( baseline_file )
                followup_files.append( followup_file )
            
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
    studies = ['ALL','ADNI1', 'ADNI2', 'ADNIGO']
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

