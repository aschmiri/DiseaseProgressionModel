#! /usr/bin/env python
# print __doc__
import argparse
import os.path
import adni_tools as adni
import csv
import sqlite3


################################################################################
#
# create_adnimerge_db
#
################################################################################
def create_adnimerge_db():
    import rpy2.robjects as robjects

    # Lead adnimerge
    robjects.r['load'](os.path.join(adni.merge_folder, 'data/adnimerge.rdata'))
    adnimerge = robjects.r['adnimerge']

    # Open database conenction
    con = sqlite3.connect(os.path.join(adni.project_folder, 'lists', 'adni.db'))
    cur = con.cursor()

    # Create table
    cur.execute('DROP TABLE IF EXISTS adnimerge_cog')
    cur.execute('CREATE TABLE adnimerge_cog (\
                                          rid INTEGER, \
                                          viscode TEXT, \
                                          apoe4 INTEGER, \
                                          cdrsb REAL, \
                                          adas11 REAL, \
                                          adas13 REAL, \
                                          faq INTEGER, \
                                          moca INTEGER, \
                                          PRIMARY KEY(rid, viscode))')

    # Print each row
    for i in range(adnimerge.nrow):
            rid = adnimerge[adnimerge.colnames.index('RID')][i]
            viscode = adnimerge[adnimerge.colnames.index('VISCODE')][i]
            apoe4 = adnimerge[adnimerge.colnames.index('APOE4')][i]
            cdrsb = adnimerge[adnimerge.colnames.index('CDRSB')][i]
            adas11 = adnimerge[adnimerge.colnames.index('ADAS11')][i]
            adas13 = adnimerge[adnimerge.colnames.index('ADAS13')][i]
            faq = adnimerge[adnimerge.colnames.index('FAQ')][i]
            moca = adnimerge[adnimerge.colnames.index('MOCA')][i]
            try:
                cur.execute("INSERT INTO adnimerge_cog VALUES (?,?,?,?,?,?,?,?)",
                            (rid, viscode, apoe4, cdrsb, adas11, adas13, faq, moca))
            except:
                print 'WARNING: ignoring double entry for', rid, viscode

    cur.execute("UPDATE adnimerge_cog SET apoe4=NULL WHERE apoe4='NA'")
    cur.execute("UPDATE adnimerge_cog SET cdrsb=NULL WHERE cdrsb='NA'")
    cur.execute("UPDATE adnimerge_cog SET adas11=NULL WHERE adas11='NA'")
    cur.execute("UPDATE adnimerge_cog SET adas13=NULL WHERE adas13='NA'")
    cur.execute("UPDATE adnimerge_cog SET faq=NULL WHERE faq='NA'")
    cur.execute("UPDATE adnimerge_cog SET moca=NULL WHERE moca='NA'")

    # Commit and close database connection
    con.commit()
    con.close()


################################################################################
#
# create_adniquery_db
#
################################################################################
def create_adniquery_db():
    # Open database conenction
    con = sqlite3.connect(os.path.join(adni.project_folder, 'lists', 'adni.db'))
    cur = con.cursor()

    # Create table
    cur.execute('DROP TABLE IF EXISTS adnimerge')
    cur.execute('CREATE TABLE adnimerge (iid INTEGER PRIMARY KEY, \
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
                                          mmse INTEGER)')

    # Read database from query file
    query_file = os.path.join(adni.project_folder, 'lists', 'query_ADNI.csv')
    with open(query_file, 'rb') as csvfile:
        entries = csv.DictReader(csvfile)
        for entry in entries:
            iid = entry['ImageUID']
            rid = entry['RID']
            viscode = entry['VISCODE']
            study_bl = entry['ORIGPROT']
            study = entry['COLPROT']
            filename = os.path.basename(entry['Files'])
            dx_bl = entry['DX.bl']
            dx = entry['DX.scan']
            scandate = entry['ScanDate']
            fs = entry['MagStrength']
            age = entry['AGE.scan']
            mmse = entry['MMSE']

            if (study == 'ADNI1' and float(fs) < 2.0) or (study != 'ADNI1' and float(fs) > 2.0):
                cur.execute("INSERT INTO adnimerge VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                            (iid, rid, viscode, study_bl, study, filename, dx_bl, dx, scandate, fs, age, mmse))

    # Commit and close database connection
    con.commit()
    con.close()


################################################################################
#
# create_years_after_diag_db
#
################################################################################
def create_years_after_diag_db():
    # Open database conenction
    con = sqlite3.connect(os.path.join(adni.project_folder, 'lists', 'adni.db'))
    cur = con.cursor()

    #
    # Read dates database
    #
    cur.execute('DROP TABLE IF EXISTS diagdates')
    cur.execute("CREATE TABLE diagdates (rid INT PRIMARY KEY, \
                                          viscode TEXT, \
                                          year_symptom INT, \
                                          year_mci INT, \
                                          year_ad INT)")

    file_dates = os.path.join(adni.project_folder, 'lists/ADNI_diag_dates.csv')

    with open(file_dates, 'rb') as csvfile:
        entries = csv.DictReader(csvfile)
        for entry in entries:
            rid = entry['Subject ID']
            viscode = entry['Visit ID']
            year_symptom = entry['Year cognitive symptoms began']
            year_mci = entry['Year MCI diagnosed']
            year_ad = entry['Year AD diagnosed']
            if viscode == 'sc' or viscode == 'v01':
                cur.execute("INSERT INTO diagdates VALUES (?,?,?,?,?)",
                            (rid, viscode, year_symptom, year_mci, year_ad))

    cur.execute("UPDATE diagdates SET year_symptom=NULL WHERE year_symptom=''")
    cur.execute("UPDATE diagdates SET year_mci=NULL WHERE year_mci=''")
    cur.execute("UPDATE diagdates SET year_ad=NULL WHERE year_ad=''")

    # Commit and close database connection
    con.commit()
    con.close()


################################################################################
#
# print_adnimerge_data_for_rid
#
################################################################################
def print_adnimerge_data_for_rid(rid=43):
    import rpy2.robjects as robjects

    # Lead database
    robjects.r['load'](os.path.join(adni.merge_folder, 'data/adnimerge.rdata'))
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
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--rid', dest='rid', type=int, default=4097)
    parser.add_argument('-q', '--create_query_db', action='store_true', default=False)
    parser.add_argument('-m', '--create_merge_db', action='store_true', default=False)
    parser.add_argument('-y', '--create_years_db', action='store_true', default=False)
    args = parser.parse_args()

    if args.create_query_db:
        create_adniquery_db()
    if args.create_merge_db:
        create_adnimerge_db()
    if args.create_years_db:
        create_years_after_diag_db()
    print_adnimerge_data_for_rid(args.rid)

    # Test DB
    con = sqlite3.connect(os.path.join(adni.project_folder, 'lists', 'adni.db'))
    cur = con.cursor()
    cur.execute("SELECT iid, rid, viscode, study, study_bl, fieldstrength, mmse, faq, moca, cdrsb FROM adnimerge JOIN adnimerge_cog USING (rid, viscode) WHERE rid = " + str(args.rid))

    rows = cur.fetchall()
    for row in rows:
        print row
