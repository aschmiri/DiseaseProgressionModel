#! /usr/bin/env python
# print __doc__

import os.path
import adni_tools as adni
import csv
import sqlite3

################################################################################
#
# create_adni_db
#
################################################################################
def create_adni_db():
    # Open database conenction
    con = sqlite3.connect( os.path.join( adni.project_folder, 'lists', 'adni.db' ) )
    cur = con.cursor()
    
    # Create table
    cur.execute( 'DROP TABLE IF EXISTS adnimerge' )
    cur.execute( 'CREATE TABLE adnimerge (iid INTEGER PRIMARY KEY, \
                                          rid INTEGER, \
                                          viscode TEXT, \
                                          study_bl TEXT, \
                                          study TEXT, \
                                          filename TEXT, \
                                          diagnosis_bl TEXT, \
                                          diagnosis TEXT, \
                                          scandate TEXT, \
                                          fieldstrength TEXT, \
                                          age REAL, \
                                          apoe4 INTEGER, \
                                          cdrsb REAL, \
                                          adas11 REAL, \
                                          adas13 REAL, \
                                          faq INTEGER, \
                                          mmse INTEGER, \
                                          moca INTEGER)' )
    
    # Read database from query file 
    query_file = os.path.join( adni.project_folder, 'lists', 'query_ADNI.csv' )
    with open( query_file, 'rb' ) as csvfile:
        entries = csv.DictReader( csvfile )
        for entry in entries:
            iid      = entry['ImageUID']
            rid      = entry['RID']
            viscode      = entry['VISCODE']
            study_bl = entry['ORIGPROT']
            study    = entry['COLPROT']
            filename = os.path.basename( entry['Files'] )
            dx_bl    = entry['DX.bl']
            dx       = entry['DX.scan']
            scandate = entry['ScanDate']
            fs       = entry['MagStrength']
            age      = entry['AGE.scan']
            apoe4    = 0#entry['APOE4']
            cdrsb    = 0#entry['CDRSB']
            adas11   = 0#entry['ADAS11']
            adas13   = 0#entry['ADAS13']
            faq      = 0#entry['FAQ']
            mmse     = entry['MMSE']
            moca     = 0#entry['MOCA']
            if (study == 'ADNI1' and float(fs) < 2.0) or (study != 'ADNI1' and float(fs) > 2.0):
                cur.execute( "INSERT INTO adnimerge VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                             (iid, rid, viscode, study_bl, study, filename, dx_bl, dx, scandate, fs, age, 
                              apoe4, cdrsb, adas11, adas13, faq, mmse, moca) )
            
    # Commit and close database connection
    con.commit()
    con.close()
    

################################################################################
#
# get_years_after_symptoms
#
################################################################################
def get_years_after_symptoms():
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    
    #
    # Read dates database 
    #
    cur.execute( "CREATE TABLE dates (rid INT PRIMARY KEY, \
                                      viscode TEXT, \
                                      year_symptom INT, \
                                      year_mci INT, \
                                      year_ad INT)" )

    file_dates = os.path.join( adni.project_folder, 'lists/ADNI_diag_dates.csv' )
    
    with open( file_dates, 'rb' ) as csvfile:
        entries = csv.DictReader( csvfile )
        for entry in entries:
            rid          = entry['Subject ID']
            viscode      = entry['Visit ID']
            year_symptom = entry['Year cognitive symptoms began']
            year_mci     = entry['Year MCI diagnosed']
            year_ad      = entry['Year AD diagnosed']
            if viscode == 'sc' or viscode == 'v01':
                cur.execute( "INSERT INTO dates VALUES (?,?,?,?,?)", 
                             (rid, viscode, year_symptom, year_mci, year_ad) )
                
    cur.execute( "UPDATE dates SET year_symptom=NULL WHERE year_symptom=''" )
    cur.execute( "UPDATE dates SET year_mci=NULL WHERE year_mci=''" )
    cur.execute( "UPDATE dates SET year_ad=NULL WHERE year_ad=''" )
        
    #
    # Read dpi database 
    #
    cur.execute( "CREATE TABLE dpis (rid INT PRIMARY KEY, \
                                     dpi FLOAT)")

    file_dates = os.path.join( adni.project_folder, 'atlas/model_1/data_sym_m24_AD.csv' )
    
    with open( file_dates, 'rb' ) as csvfile:
        entries = csv.DictReader( csvfile )
        for entry in entries:
            rid      = entry['RID']
            dpi      = entry['DPI']
            cur.execute( "INSERT INTO dpis VALUES (?,?)", 
                         (rid, dpi) )
    
    #
    # Read summary database 
    #
    cur.execute( "CREATE TABLE summary (rid INT PRIMARY KEY, \
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
                                        moca INT)" )

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
                    cur.execute( "INSERT INTO summary VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                                 (rid, filename, dx, scandate, age, apoe4, cdrsb, adas11, adas13, faq, mmse, moca) )

    add_summary( '/vol/medic01/users/aschmidt/projects/Data/ADNI/data/ADNI1/native/images', 
                 os.path.join( adni.project_folder, 'lists/query_ADNI1.csv' ) )
    add_summary( '/vol/medic01/users/aschmidt/projects/Data/ADNI/data/ADNI2/native/images', 
                 os.path.join( adni.project_folder, 'lists/query_ADNI2.csv' ) )
    add_summary( '/vol/medic01/users/aschmidt/projects/Data/ADNI/data/ADNIGO/native/images', 
                 os.path.join( adni.project_folder, 'lists/query_ADNIGO.csv' ) )
    
    cur.execute( "UPDATE summary SET apoe4=NULL WHERE apoe4='NA'" )
    cur.execute( "UPDATE summary SET cdrsb=NULL WHERE cdrsb='NA'" )
    cur.execute( "UPDATE summary SET adas11=NULL WHERE adas11='NA'" )
    cur.execute( "UPDATE summary SET adas13=NULL WHERE adas13='NA'" )
    cur.execute( "UPDATE summary SET faq=NULL WHERE faq='NA'" )
    cur.execute( "UPDATE summary SET mmse=NULL WHERE mmse='NA'" )
    cur.execute( "UPDATE summary SET moca=NULL WHERE moca='NA'" )
    
    cur.execute( "SELECT rid, scandate, year_symptom, dpi, cdrsb, adas11, adas13, faq, mmse, moca \
                  FROM summary NATURAL JOIN dates NATURAL JOIN dpis" )
    
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
    robjects.r['load']( os.path.join( adni.merge_folder, 'data/adnimerge.rdata' ) )
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
    #get_years_after_symptoms()
    print_data_for_rid( 54 )
    #create_adni_db()
    
    # Test DB
    con = sqlite3.connect( os.path.join( adni.project_folder, 'lists', 'adni.db' ) )
    cur = con.cursor()
    cur.execute( "SELECT iid, rid, viscode, study, study_bl, fieldstrength, mmse FROM adnimerge WHERE rid = 54" )
    
    rows = cur.fetchall()
    for row in rows:
        print row

