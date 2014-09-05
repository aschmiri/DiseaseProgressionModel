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
from vgam.synthmodel import SynthModel


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
        parser.add_argument('method', choices=['cog', 'reg', 'long', 'cons', 'graph', 'mbl', 'synth', 'all'])
        parser.add_argument('-i', '--iteration', type=int, default=0, help='the refinement iteration')
        parser.add_argument('-b', '--biomarker_name', default=None, help='name of the biomarker to be plotted')
        parser.add_argument('--trans', type=str, default='sym', help='the transformation model, e.g. ffd, svffd, sym, or ic (regbased only)')
        parser.add_argument('--spacing', type=str, default='5', help='the transformation spacing (regbased only)')
        return parser

    ############################################################################
    #
    # get_data_handler()
    #
    ############################################################################
    @staticmethod
    def get_data_handler(args=None, iteration=None):
        if args is not None and args.method == 'synth':
            return SynthDataHandler(args=args, iteration=iteration)
        else:
            return DataHandler(args=args, iteration=iteration)

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
        if args is None:
            self._single_biomarker = None
        else:
            self._single_biomarker = args.biomarker_name
        self._biomarker_sets = {'cog': adni.cog_score_names,
                                'reg': adni.volume_names,
                                'long': adni.volume_names,
                                'cons': adni.volume_names,
                                'graph': adni.volume_names_essential,
                                'mbl': adni.manifold_coordinate_names,
                                'all': adni.biomarker_names_essential}

        # Set data folders
        self._model_folders = {'cog': os.path.join(adni.project_folder, 'models', 'cog'),
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

        self._data_files = {'meta': os.path.join(adni.project_folder, 'lists/metadata.csv'),
                            'cog': os.path.join(adni.project_folder, 'lists/metadata.csv'),
                            'mbl': os.path.join(adni.project_folder, 'lists/manifold_features.csv'),
                            'reg': os.path.join(adni.project_folder, 'lists/volumes_segbased_' + trans + '_' + spacing + 'mm.csv'),
                            'long': os.path.join(adni.project_folder, 'lists/volumes_segbased_longitudinal.csv'),
                            'cons': os.path.join(adni.project_folder, 'lists/volumes_segbased_consistent.csv'),
                            'graph': os.path.join(adni.project_folder, 'lists/volumes_segbased_graphcut.csv')}

        # Set method
        if args is None:
            self._method = 'all'
            self._volume_method = 'graph'
        else:
            self._method = args.method
            if args.method == 'reg':
                self._volume_method = 'reg'
            elif args.method == 'long':
                self._volume_method = 'long'
            elif args.method == 'cons':
                self._volume_method = 'cons'
            else:
                self._volume_method = 'graph'

        # Set iteration
        if iteration is not None:
            self._iteration = iteration
        elif args is not None:
            self._iteration = args.iteration
        else:
            self._iteration = 0

    ############################################################################
    #
    # get_biomarker_set()
    #
    ############################################################################
    def get_biomarker_set(self, method=None):
        ''' Get the right set of biomarkers for the given method.'''
        if self._single_biomarker is not None:
            return [self._single_biomarker]
        method = self._method if method is None else method
        return self._biomarker_sets[method]

    ############################################################################
    #
    # get_data_file()
    #
    ############################################################################
    def get_data_file(self, method=None):
        ''' Get the right data file for the given method. '''
        method = self._method if method is None else method
        return self._data_files[method]

    ############################################################################
    #
    # get_data_file_for_biomarker()
    #
    ############################################################################
    def get_data_file_for_biomarker(self, biomarker):
        ''' Get the right data file for the given biomarker. '''
        return self._data_files[self._get_method_for_biomarker(biomarker)]

    ############################################################################
    #
    # get_model_folders()
    #
    ############################################################################
    def get_model_folder(self, method=None, previous_iteration=False):
        ''' Get the right data file for the given method and iteration. '''
        method = self._method if method is None else method
        iteration = self._iteration - 1 if previous_iteration else self._iteration
        iteration_folder = 'it_{0}'.format(iteration)
        return adni.make_dir(self._model_folders[method], iteration_folder)

    ############################################################################
    #
    # get_model_folder_for_biomarker()
    #
    ############################################################################
    def get_model_folder_for_biomarker(self, biomarker, iteration=None):
        ''' Get the right data folders for the given biomarker. '''
        iteration = self._iteration if iteration is None else iteration
        base_folder = self._model_folders[self._get_method_for_biomarker(biomarker)]
        iteration_folder = 'it_{0}'.format(iteration)
        return adni.make_dir(base_folder, iteration_folder)

    ############################################################################
    #
    # get_model_file()
    #
    ############################################################################
    def get_model_file(self, biomarker, iteration=None):
        ''' Get the right model file for the given biomarker. '''
        iteration = self._iteration if iteration is None else iteration
        model_folder = self.get_model_folder_for_biomarker(biomarker, iteration=iteration)
        return os.path.join(model_folder, biomarker.replace(' ', '_') + '_model.csv')

    ############################################################################
    #
    # get_samples_file()
    #
    ############################################################################
    def get_samples_file(self, biomarker, iteration=None):
        ''' Get the right model file for the given biomarker. '''
        iteration = self._iteration if iteration is None else iteration
        model_folder = self.get_model_folder_for_biomarker(biomarker, iteration=iteration)
        return os.path.join(model_folder, biomarker.replace(' ', '_') + '_samples.csv')

    ############################################################################
    #
    # _get_method_for_biomarker()
    #
    ############################################################################
    def _get_method_for_biomarker(self, biomarker):
        if biomarker in adni.cog_score_names:
            return 'cog'

        elif biomarker in adni.volume_names:
            return self._volume_method

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
        measurements = self._get_metadata_as_dict(select_converters=select_converters)
        measurements = self.update_measurements_with_biomarker_values(measurements, biomarkers=biomarkers, no_regression=no_regression)
        if complete:
            measurements = self._ensure_complete_measurements(measurements, biomarkers=biomarkers)

        # Return measurements
        return measurements

    ############################################################################
    #
    # _ensure_complete_measurements()
    #
    ############################################################################
    def _ensure_complete_measurements(self, measurements, biomarkers=None):
        biomarkers = self.get_biomarker_set() if biomarkers is None else biomarkers

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
    # _get_metadata_as_dict()
    #
    ############################################################################
    def _get_metadata_as_dict(self, select_converters=True):
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
        biomarkers = self.get_biomarker_set() if biomarkers is None else biomarkers
        biomarker_values = self._get_biomarker_values_as_dict(measurements, biomarkers=biomarkers, no_regression=no_regression)

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
    # _get_biomarker_values_as_dict()
    #
    ############################################################################
    def _get_biomarker_values_as_dict(self, metadata, biomarkers=None, no_regression=False):
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
        biomarkers = self.get_biomarker_set() if biomarkers is None else biomarkers

        NO_REGRESSION_STR = '{0} no regression'
        print log.INFO, 'Collecting biomarker values...'
        values = {}

        for biomarker in biomarkers:
            data_file = self.get_data_file(method=self._get_method_for_biomarker(biomarker))
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
                    regressed_vals = self._age_regression(ages, vals)

                    for rid, viscode, value in zip(rids, viscodes, regressed_vals):
                        values[rid][viscode].update({biomarker: value})

        return values

    ############################################################################
    #
    # _age_regression()
    #
    ############################################################################
    def _age_regression(self, timepoints, values):
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
    # get_mmse_decline_as_dict()
    #
    ############################################################################
    def get_mmse_decline_as_dict(self):
        ''' Return all a dictionary that indicates for each RID the MMSE decline
        of the subject between baseline and m24 visits.

        Arguments:
        measurements -- the measurements as a dictionary

        Returns:
        A dictionary with the following structure:
        { <rid> : <decline>
          ... }
        '''
        measurements = self._get_metadata_as_dict(select_converters=False)
        measurements = self.update_measurements_with_biomarker_values(measurements, biomarkers=['MMSE'])

        declines = {}
        for rid in measurements:
            try:
                mmse_bl = measurements[rid]['bl']['MMSE']
                mmse_24 = measurements[rid]['m24']['MMSE']
                decline = mmse_bl - mmse_24
                declines.update({rid: decline})
            except Exception:
                # TODO: debug message
                pass

        return declines

    ############################################################################
    #
    # get_rcd_as_dict()
    #
    ############################################################################
    def get_rcd_as_dict(self):
        ''' Return all a dictionary that indicates for each RID if the subject
        is classified as RCD (rapid cognitive decline).

        Arguments:
        measurements -- the measurements as a dictionary

        Returns:
        A dictionary with the following structure:
        { <rid> : [True|False]
          ... }
        '''
        measurements = self._get_metadata_as_dict(select_converters=False)
        measurements = self.update_measurements_with_biomarker_values(measurements, biomarkers=['MMSE'])

        rcds = {}
        for rid in measurements:
            try:
                mmse_bl = measurements[rid]['bl']['MMSE']
                mmse_24 = measurements[rid]['m24']['MMSE']
                rcd = True if (mmse_24 - mmse_bl) < -7 else False
                rcds.update({rid: rcd})
            except Exception:
                # TODO: debug message
                pass

        return rcds


class SynthDataHandler(DataHandler):
    '''
    TODO: classdocs
    '''

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
        args.iteration -- the iteration of the fitting
        '''
        # Set biomarker sets
        self._single_biomarker = None
        self._biomarker_sets = {'synth': SynthModel.get_biomarker_names()}

        # Set data folders
        self._model_folders = {'synth': os.path.join(adni.project_folder, 'models', 'synth')}

        # Set data files
        self._data_files = {'meta': os.path.join(adni.project_folder, 'lists/synthdata.csv'),
                            'synth': os.path.join(adni.project_folder, 'lists/synthdata.csv')}

        # Set method
        self._method = 'synth'
        self._volume_method = 'synth'

        # Set iteration
        if iteration is not None:
            self._iteration = iteration
        elif args is not None:
            self._iteration = args.iteration
        else:
            self._iteration = 0

    ############################################################################
    #
    # _get_metadata_as_dict()
    #
    ############################################################################
    def _get_metadata_as_dict(self, select_converters=True):
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

                # Get progress
                progress = row['progress']
                metadata[rid][viscode].update({'progress': progress})

        return metadata

    ############################################################################
    #
    # _get_biomarker_values_as_dict()
    #
    ############################################################################
    def _get_biomarker_values_as_dict(self, metadata, biomarkers=None, no_regression=False):
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
        biomarkers = self.get_biomarker_set() if biomarkers is None else biomarkers

        print log.INFO, 'Collecting biomarker values...'
        values = {}

        for biomarker in biomarkers:
            data_file = self.get_data_file(method=self._get_method_for_biomarker(biomarker))
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
                        values[rid][viscode].update({biomarker: value})

        return values

    ############################################################################
    #
    # _get_method_for_biomarker()
    #
    ############################################################################
    def _get_method_for_biomarker(self, biomarker):
        return 'synth'

    ############################################################################
    #
    # get_mmse_decline_as_dict()
    #
    ############################################################################
    def get_mmse_decline_as_dict(self):
        print log.ERROR, 'MMSE decline not available for synthetic data!'
        return None

    ############################################################################
    #
    # get_rcd_as_dict()
    #
    ############################################################################
    def get_rcd_as_dict(self):
        print log.ERROR, 'RCD not available for synthetic data!'
        return None
