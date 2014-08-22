'''
A class to provide progression modelling functionality.

@author:     Alexander Schmidt-Richberg
@copyright:  2014 Imperial College London. All rights reserved.
@contact:    a.schmidt-richberg@imperial.ac.uk
'''
import os.path
import datetime
import csv
import numpy as np
from common import log as log
from common import adni_tools as adni


class DataHandler(object):
    '''
    TODO: classdocs
    '''
    ############################################################################
    #
    # add_arguments()
    #
    ############################################################################
    @staticmethod
    def add_arguments(parser):
        parser.add_argument('method', choices=['cog', 'reg', 'long', 'cons', 'graph', 'mbl'])
        parser.add_argument('-i', '--iteration', type=int, default=0, help='the refinement iteration')
        parser.add_argument('-b', '--biomarker_name', default=None, help='name of the biomarker to be plotted')
        parser.add_argument('--trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic (regbased only)')
        parser.add_argument('--spacing', type=str, default='5', help='the transformation spacing (regbased only)')
        return parser

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, args=None, iteration=None):
        '''
        Initialise the right data files for the given arguments.

        Arguments:
        args -- command line arguments with:
        args.method -- the method, choice of ['cog', 'reg', 'long', 'cons', 'graph', 'mbl']
        args.iteration -- the iteration of the fitting
        args.trans -- the transformation if method == 'reg'
        args.spacing -- the spacing if method == 'reg'
        '''
        # Set biomarker sets
        self.__biomarker_sets = {'cog': adni.cog_score_names,
                                 'reg': adni.volume_names,
                                 'long': adni.volume_names,
                                 'cons': adni.volume_names,
                                 'graph': adni.volume_names_essential,
                                 'mbl': adni.manifold_coordinate_names,
                                 'all': adni.biomarker_names}

        # Set data folders
        self.__model_folders = {'cog': os.path.join(adni.project_folder, 'models', 'cog'),
                                'mbl': os.path.join(adni.project_folder, 'models', 'mbl'),
                                'reg': os.path.join(adni.project_folder, 'models', 'reg'),
                                'long': os.path.join(adni.project_folder, 'models', 'long'),
                                'cons': os.path.join(adni.project_folder, 'models', 'cons'),
                                'graph': os.path.join(adni.project_folder, 'models', 'graph')}

        # Set data files
        if args is None:
            trans = 'sym'
            spacing = '5'
        else:
            trans = args.trans
            spacing = args.spacing

        self.__data_files = {'meta': os.path.join(adni.project_folder, 'lists/metadata.csv'),
                             'cog': os.path.join(adni.project_folder, 'lists/metadata.csv'),
                             'mbl': os.path.join(adni.project_folder, 'lists/manifold_features.csv'),
                             'reg': os.path.join(adni.project_folder, 'lists/volumes_segbased_' + trans + '_' + spacing + 'mm.csv'),
                             'long': os.path.join(adni.project_folder, 'lists/volumes_segbased_longitudinal.csv'),
                             'cons': os.path.join(adni.project_folder, 'lists/volumes_segbased_consistent.csv'),
                             'graph': os.path.join(adni.project_folder, 'lists/volumes_segbased_graphcut.csv')}

        # Set method
        if args is None:
            self.__method = 'all'
            self.__volume_method = 'graph'
        else:
            self.__method = args.method
            if args.method == 'reg':
                self.__volume_method = 'reg'
            elif args.method == 'long':
                self.__volume_method = 'long'
            elif args.method == 'cons':
                self.__volume_method = 'cons'
            else:
                self.__volume_method = 'graph'

        # Set iteration
        if iteration is not None:
            self.__iteration = iteration
        elif args is not None:
            self.__iteration = args.iteration
        else:
            self.__iteration = 0

    ############################################################################
    #
    # get_biomarker_set()
    #
    ############################################################################
    def get_biomarker_set(self, method=None):
        ''' Get the right set of biomarkers for the given method.'''
        if method is None:
            method = self.__method
        return self.__biomarker_sets[method]

    ############################################################################
    #
    # get_data_file()
    #
    ############################################################################
    def get_data_file(self, method=None):
        ''' Get the right data file for the given method. '''
        if method is None:
            method = self.__method
        return self.__data_files[method]

    ############################################################################
    #
    # get_data_file_for_biomarker()
    #
    ############################################################################
    def get_data_file_for_biomarker(self, biomarker):
        ''' Get the right data file for the given biomarker. '''
        return self.__data_files[self.__get_method_for_biomarker(biomarker)]

    ############################################################################
    #
    # get_model_folders()
    #
    ############################################################################
    def get_model_folder(self, method=None, previous_iteration=False):
        ''' Get the right data file for the given method and iteration. '''
        if method is None:
            method = self.__method

        if previous_iteration:
            iteration = self.__iteration - 1
        else:
            iteration = self.__iteration

        iteration_folder = 'it_{0}'.format(iteration)
        return adni.make_dir(self.__model_folders[method], iteration_folder)

    ############################################################################
    #
    # get_model_folder_for_biomarker()
    #
    ############################################################################
    def get_model_folder_for_biomarker(self, biomarker, iteration=None):
        ''' Get the right data folders for the given biomarker. '''
        if iteration is None:
            iteration = self.__iteration

        base_folder = self.__model_folders[self.__get_method_for_biomarker(biomarker)]

        iteration_folder = 'it_{0}'.format(iteration)
        return adni.make_dir(base_folder, iteration_folder)

    ############################################################################
    #
    # get_model_file()
    #
    ############################################################################
    def get_model_file(self, biomarker, iteration=None):
        ''' Get the right model file for the given biomarker. '''
        if iteration is None:
            iteration = self.__iteration
        model_folder = self.get_model_folder_for_biomarker(biomarker, iteration=iteration)
        return os.path.join(model_folder, biomarker.replace(' ', '_') + '_model.csv')

    ############################################################################
    #
    # get_samples_file()
    #
    ############################################################################
    def get_samples_file(self, biomarker, iteration=None):
        ''' Get the right model file for the given biomarker. '''
        if iteration is None:
            iteration = self.__iteration
        model_folder = self.get_model_folder_for_biomarker(biomarker, iteration=iteration)
        return os.path.join(model_folder, biomarker.replace(' ', '_') + '_samples.csv')

    ############################################################################
    #
    # __get_method_for_biomarker()
    #
    ############################################################################
    def __get_method_for_biomarker(self, biomarker):
        if biomarker in adni.cog_score_names:
            return 'cog'

        elif biomarker in adni.volume_names:
            return self.__volume_method

        elif biomarker in adni.manifold_coordinate_names:
            return 'mbl'

        else:
            print log.ERROR, 'Tag cannot be determined.'
            return None

    ############################################################################
    #
    # get_measurements_as_dict()
    #
    ############################################################################
    def get_measurements_as_dict(self, biomarkers=None, select_converters=True, no_regression=False, complete=False):
        ''' Return all subjects measurements as a dictionary.

        Arguments:
        select_converters -- only select MCI -> AD converters
        no_regression -- do not perform age regression

        Returns:
        A dictionary with the following structure:
        { <rid> : { <viscode> : { DX.scan : <diagnosis> }
                                { AGE.scan : <age in years> }
                                { scandate : <date of scan> }
                                { scantime : <days after bl> }
                                { progress : <days relative to conversion> }
                                { <biomarker1> : <volume> }
                                ... }
                  { <viscode> : ... }
          <rid> : ... }
        '''
        # Read data from lists
        measurements = self.__get_metadata_as_dict(select_converters=select_converters)
        measurements = self.update_measurements_with_biomarker_values(measurements, biomarkers=biomarkers, no_regression=no_regression)
        if complete:
            measurements = self.__ensure_complete_measurements(measurements, biomarkers=biomarkers)

        # Return measurements
        return measurements

    ############################################################################
    #
    # __ensure_complete_measurements()
    #
    ############################################################################
    def __ensure_complete_measurements(self, measurements, biomarkers=None):
        if biomarkers is None:
            biomarkers = self.get_biomarker_set()

        for rid, visit in measurements.items():
            drop_rid = False
            for _, visdata in visit.items():
                for biomarker in biomarkers:
                    if biomarker not in visdata:
                        drop_rid = True
                        break
            if drop_rid:
                measurements.pop(rid)

        return measurements

    ############################################################################
    #
    # __get_metadata_as_dict()
    #
    ############################################################################
    def __get_metadata_as_dict(self, select_converters=True):
        ''' Return all subjects metadata as a dictionary.

        Arguments:
        data_file -- the names of the csv file containing the data
        select_converters -- only select MCI -> AD converters

        Returns:
        A dictionary with the following structure:
        { <rid> : { <viscode> : { DX.scan : <diagnosis> }
                                { AGE.scan : <age in years> }
                                { scandate : <date of scan> }
                                { scantime : <days after bl> }
                                { progress : <days relative to conversion> } }
                  { <viscode> : ... }
          <rid> : ... }
        '''
        print log.INFO, 'Collecting metadata...'
        metadata = {}

        data_file = self.get_data_file(method='meta')
        if not os.path.isfile(data_file):
            print log.ERROR, 'Data file dies not exist:', data_file
            return metadata

        #
        # Get all measurements from CSV file
        with open(data_file, 'rb') as csvfile:
            rows = csv.DictReader(csvfile)
            for row in rows:
                # Get rid
                rid = int(row['RID'])
                if rid not in metadata:
                    metadata.update({rid: {}})

                # Get scan time
                viscode = row['VISCODE']
                if viscode in metadata[rid]:
                    print log.WARNING, 'Entry already exists {0} ({1}). Skipping.'.format(rid, viscode)
                    break
                metadata[rid].update({viscode: {}})

                # Get scan date
                scandate = datetime.datetime.strptime(row['ScanDate'], "%Y-%m-%d").date()
                metadata[rid][viscode].update({'scandate': scandate})

                # Get age
                metadata[rid][viscode].update({'AGE.scan': float(row['AGE.scan'])})

                # Get factor
                metadata[rid][viscode].update({'FactorMNI': float(row['FactorMNI'])})

                # Get diagnosis as numerical value
                dx_str = row['DX.scan']
                if dx_str == 'AD':
                    dx = 1.0
                elif dx_str == 'MCI':
                    dx = 0.5
                elif dx_str == 'CN':
                    dx = 0.0
                else:
                    print log.ERROR, 'Invalid diagnosis:', viscode
                    break
                metadata[rid][viscode].update({'DX.scan': dx})

        #
        # Add time relative to point of conversion to each data set
        if select_converters:
            print log.INFO, 'Selecting converters...'
            valid_rids = []
            for rid, rid_data in metadata.items():
                data_viscode = []
                data_scantime = []
                data_diagnosis = []

                if 'bl' not in metadata[rid]:
                    print log.WARNING, 'No bl scan for subject {0}!'.format(rid)
                    continue

                bl_date = metadata[rid]['bl']['scandate']
                for viscode, scan_data in rid_data.items():
                    fu_date = metadata[rid][viscode]['scandate']
                    scantime = (fu_date - bl_date).days

                    data_viscode.append(viscode)
                    data_scantime.append(scantime)
                    data_diagnosis.append(scan_data['DX.scan'])

                data_viscode = np.array(data_viscode)
                data_scantime = np.array(data_scantime)
                data_diagnosis = np.array(data_diagnosis)

                args = np.argsort(data_scantime)
                data_viscode = data_viscode[args]
                data_scantime = data_scantime[args]
                data_diagnosis = data_diagnosis[args]

                if data_diagnosis[-1] == 1.0 and data_diagnosis[0] == 0.5:
                    valid_rids.append(rid)
                    scantime_prev = data_scantime[0]
                    for diagnosis, scantime in zip(data_diagnosis, data_scantime):
                        if diagnosis == 1.0:
                            time_convert = scantime_prev + (scantime - scantime_prev) / 2
                            break
                        else:
                            scantime_prev = scantime

                    for viscode, scantime in zip(data_viscode, data_scantime):
                        metadata[rid][viscode].update({'scantime': scantime})
                        metadata[rid][viscode].update({'progress': scantime - time_convert})

            #
            # Remove non-valid rids
            metadata = {rid: value for rid, value in metadata.items() if rid in valid_rids}

        return metadata

    ############################################################################
    #
    # update_measurements_with_biomarker_values()
    #
    ############################################################################
    def update_measurements_with_biomarker_values(self, measurements, biomarkers=None, no_regression=False):
        ''' Update the metadata dictionary with the biomarker values.

        Arguments:
        metadata -- the metadata
        data_file -- the data file containing the biomarker values
        biomarker_names -- the biomarkers to be read
        no_regression -- do not perform age regression
        '''
        if biomarkers is None:
            biomarkers = self.get_biomarker_set()

        biomarker_values = self.__get_biomarker_values_as_dict(measurements, biomarkers=biomarkers, no_regression=no_regression)

        # Update metadata with feature values
        for rid in measurements:
            for viscode in measurements[rid]:
                for biomarker in biomarkers:
                    try:
                        value = biomarker_values[rid][viscode][biomarker]
                        measurements[rid][viscode].update({biomarker: value})
                    except:
                        # TODO: Debug message
                        pass

        return measurements

    ############################################################################
    #
    # __get_biomarker_values_as_dict()
    #
    ############################################################################
    def __get_biomarker_values_as_dict(self, metadata, biomarkers=None, no_regression=False):
        ''' Return all measurements as a dictionary.

        Arguments:
        data_file -- the names of the csv file containing the data
        biomarker_names -- the biomarkers to be read
        no_regression -- do not perform age regression

        Returns:
        A dictionary with the following structure:
        { <rid> : { <viscode> : { <biomarker1> : <value> }
                                { <biomarker2> : <value> } ... }
                  { <viscode> : ... }
          <rid> : ... }
        '''
        if biomarkers is None:
            biomarkers = self.get_biomarker_set()

        NO_REGRESSION_STR = '{0} no regression'
        print log.INFO, 'Collecting biomarker values...'
        values = {}

        for biomarker in biomarkers:
            data_file = self.get_data_file(method=self.__get_method_for_biomarker(biomarker))
            if data_file is None:
                print log.ERROR, 'Data file not found:', data_file
                continue

            #
            # Get all measurements from CSV file
            with open(data_file, 'rb') as csvfile:
                rows = csv.DictReader(csvfile)
                for row in rows:
                    # Get rid
                    rid = int(row['RID'])
                    viscode = row['VISCODE']
                    if rid in metadata and viscode in metadata[rid]:
                        # Initialise dictionary
                        if rid not in values:
                            values.update({rid: {}})
                        if viscode not in values[rid]:
                            values[rid].update({viscode: {}})

                        # Get value
                        value = adni.safe_cast(row[biomarker])

                        # Normalise volumes with mni transformation
                        if biomarker in adni.volume_names:
                            factor = metadata[rid][viscode]['FactorMNI']
                            value /= factor

                        if no_regression:
                            values[rid][viscode].update({biomarker: value})
                        else:
                            values[rid][viscode].update({NO_REGRESSION_STR.format(biomarker): value})

        if not no_regression:
            print log.INFO, 'Performing age regression...'
            for biomarker in biomarkers:
                rids = []
                viscodes = []
                vals = []
                ages = []
                for rid, visits in values.items():
                    for viscode in visits:
                        # Get age
                        try:
                            age = metadata[rid][viscode]['AGE.scan']
                        except:
                            age = None

                        if age is not None and NO_REGRESSION_STR.format(biomarker) in visits[viscode]:
                            value = visits[viscode][NO_REGRESSION_STR.format(biomarker)]
                            if value is not None:
                                ages.append(age)
                                vals.append(value)
                                rids.append(rid)
                                viscodes.append(viscode)

                if len(ages) > 1:
                    regressed_vals = self.__age_regression(ages, vals)

                    for rid, viscode, value in zip(rids, viscodes, regressed_vals):
                        values[rid][viscode].update({biomarker: value})

        return values

    ############################################################################
    #
    # __age_regression()
    #
    ############################################################################
    def __age_regression(self, timepoints, values):
        ''' Perform age regression. '''
        timepoints = np.array(timepoints)
        values = np.array(values)
        mean_val = np.mean(values)

        import scipy.stats as stats
        slope, intercept, _, _, _ = stats.linregress(timepoints, values)
        values = values - (timepoints * slope + intercept) + mean_val

        return values

    ############################################################################
    #
    # get_rcd_as_dict()
    #
    ############################################################################
    def get_rcd_as_dict(self, measurements):
        ''' Return all a dictionary that indicates for each RID if the subject
        is classified as RCD (rapid cognitive decline.

        Arguments:
        measurements -- the measurements as a dictionary

        Returns:
        A dictionary with the following structure:
        { <rid> : [True|False]
          ... }
        '''
        rcds = {}
        for rid in measurements:
            try:
                mmse_bl = measurements[rid]['bl']['MMSE']
                mmse_24 = measurements[rid]['m24']['MMSE']
                # rcd = True if (mmse_24 - mmse_bl) < -7 else False
                decline = mmse_bl - mmse_24
                rcds.update({rid: decline})
            except Exception:
                # TODO: debug message
                pass

        return rcds
