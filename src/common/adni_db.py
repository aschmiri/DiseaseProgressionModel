#! /usr/bin/env python
# print __doc__

import os.path
import adni_tools as adni

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
    print_data_for_rid( 1226 )

