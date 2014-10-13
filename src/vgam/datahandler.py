"""
A class to provide reading and writing data for progression modelling.

:author:     Alexander Schmidt-Richberg
:copyright:  2014 Imperial College London. All rights reserved.
:contact:    a.schmidt-richberg@imperial.ac.uk
"""
import os.path
import datetime
import csv
import numpy as np
from common import log as log
from common import adni_tools as adni
from vgam.synthmodel import SynthModel


class DataHandler(object):
    """
    TODO: classdocs
    """
    ############################################################################
    #
    # add_arguments()
    #
    ############################################################################
    @staticmethod
    def add_arguments(parser):
        parser.add_argument('-m', '--method', choices=['cog', 'long', 'cons', 'mbl', 'img', 'all', 'synth'], default='all', help='the method to collect data for')
        parser.add_argument('-i', '--iteration', type=int, default=0, help='the refinement iteration')
        parser.add_argument('-b', '--biomarkers_name', nargs='+', default=None, help='name of the biomarker to be plotted')
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
        """
        Initialise the right data files for the given arguments.

        :param args: command line arguments with:
        :param str args.method: the method, choice of ['cog', 'long', 'cons', 'mbl', 'img', 'all']
        :param int args.iteration: the iteration of the fitting
        """
        # Set biomarker sets
        if args is None:
            self._biomarker_subset = None
        else:
            self._biomarker_subset = args.biomarkers_name
        self._biomarker_sets = {'cog': adni.cog_score_names,
                                'long': adni.structure_names,
                                'cons': adni.structure_names,
                                'mbl': adni.manifold_coordinate_names,
                                'img': adni.image_biomarker_names,
                                'all': adni.biomarker_names}

        # Set data folders
        self._model_folders = {'cog': os.path.join(adni.project_folder, 'models', 'cog'),
                               'mbl': os.path.join(adni.project_folder, 'models', 'mbl'),
                               'long': os.path.join(adni.project_folder, 'models', 'long'),
                               'cons': os.path.join(adni.project_folder, 'models', 'cons')}

        # Set data files
        self._data_files = {'meta': os.path.join(adni.project_folder, 'lists/metadata.csv'),
                            'cog': os.path.join(adni.project_folder, 'lists/metadata.csv'),
                            'mbl': os.path.join(adni.project_folder, 'lists/manifold_features.csv'),
                            'long': os.path.join(adni.project_folder, 'lists/volumes_segbased_longitudinal.csv'),
                            'cons': os.path.join(adni.project_folder, 'lists/volumes_segbased_consistent.csv')}

        # Set method
        if args is None:
            self._method = 'all'
            self._volume_method = 'long'
        else:
            self._method = args.method
            if args.method == 'long':
                self._volume_method = 'long'
            elif args.method == 'cons':
                self._volume_method = 'cons'
            else:
                self._volume_method = 'long'

        # Set iteration
        if iteration is not None:
            self._iteration = iteration
        elif args is not None:
            self._iteration = args.iteration
        else:
            self._iteration = 0

        self._diagnosis_code = {'CN': 0.0, 'EMCI': 0.25, 'MCI': 0.5, 'LMCI': 0.75, 'AD': 1.0}

    ############################################################################
    #
    # get_biomarker_set()
    #
    ############################################################################
    def get_biomarker_set(self, method=None):
        """ Get the right set of biomarkers for the given method."""
        if self._biomarker_subset is not None:
            return self._biomarker_subset
        method = self._method if method is None else method
        return self._biomarker_sets[method]

    ############################################################################
    #
    # get_data_file()
    #
    ############################################################################
    def get_data_file(self, method=None):
        """ Get the right data file for the given method. """
        method = self._method if method is None else method
        return self._data_files[method]

    ############################################################################
    #
    # get_data_file_for_biomarker()
    #
    ############################################################################
    def get_data_file_for_biomarker(self, biomarker):
        """ Get the right data file for the given biomarker. """
        return self._data_files[self._get_method_for_biomarker(biomarker)]

    ############################################################################
    #
    # get_model_folders()
    #
    ############################################################################
    def get_model_folder(self, method=None, previous_iteration=False):
        """ Get the right data file for the given method and iteration. """
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
        """ Get the right data folders for the given biomarker. """
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
        """ Get the right model file for the given biomarker. """
        iteration = self._iteration if iteration is None else iteration
        model_folder = self.get_model_folder_for_biomarker(biomarker, iteration=iteration)
        return os.path.join(model_folder, biomarker.replace(' ', '_') + '_model.csv')

    ############################################################################
    #
    # get_samples_file()
    #
    ############################################################################
    def get_samples_file(self, biomarker, iteration=None):
        """ Get the right model file for the given biomarker. """
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
        elif biomarker in adni.structure_names or biomarker in adni.structure_names_complete:
            return self._volume_method
        elif biomarker in adni.manifold_coordinate_names:
            return 'mbl'
        else:
            print log.ERROR, 'Tag cannot be determined for {0}!'.format(biomarker)
            return None

    ############################################################################
    #
    # get_measurements_as_dict()
    #
    ############################################################################
    def get_measurements_as_dict(self, min_visits=0, visits=None, biomarkers=None,
                                 select_training_set=False, select_test_set=False,
                                 select_complete=False, no_regression=False):
        """ Return all subjects measurements as a dictionary.

        Arguments:
        :param int min_visits: minimal number of visits that have to be available
        :param list visits: list of specific visits
        :param list biomarkers: list of specific biomarkers
        :param bool select_training_set: only select MCI -> AD converters
        :param bool select_test_set: select the test subjects
        :param bool select_complete: select subjects with complete biomarker measurements
        :param bool no_regression: do not perform age regression

        Returns:
        :return: a dictionary with the following structure:
        { <rid> : { <viscode> : { DX.scan : <diagnosis> }
                                { AGE.scan : <age in years> }
                                { scandate : <date of scan> }
                                { scantime : <days after bl> }
                                { progress : <days relative to conversion> }
                                { <biomarker1> : <volume> }
                                ... }
                  { <viscode> : ... }
          <rid> : ... }
        :rtype: dict
        """
        # Read data from lists
        measurements = self._get_metadata_as_dict()
        measurements = self.update_measurements_with_biomarker_values(measurements,
                                                                      biomarkers=biomarkers,
                                                                      no_regression=no_regression)

        # Select specific subset of data. this has to be done after age regression!
        if select_training_set:
            measurements = self._select_training_set(measurements)
        if select_test_set:
            measurements = self._select_test_set(measurements)
        if visits is not None or min_visits > 0:
            measurements = self._select_visits(measurements, min_visits=min_visits, visits=visits)
        if select_complete:
            measurements = self._select_complete_measurements(measurements, biomarkers=biomarkers, visits=visits)

        # Return measurements=
        return measurements

    ############################################################################
    #
    # _get_metadata_as_dict()
    #
    ############################################################################
    def _get_metadata_as_dict(self):
        """ Return all subjects metadata as a dictionary.

        :param int min_visits: minimal number of visits

        :return: a dictionary with the following structure:
        { <rid> : { <viscode> : { DX.scan : <diagnosis> }
                                { AGE.scan : <age in years> }
                                { scandate : <date of scan> }
                                { scantime : <days after bl> }
                                { progress : <days relative to conversion> } }
                  { <viscode> : ... }
          <rid> : ... }
        :rtype dict
        """
        print log.INFO, 'Collecting metadata...'
        metadata = {}

        data_file = self.get_data_file(method='meta')
        if not os.path.isfile(data_file):
            print log.ERROR, 'Data file dies not exist:', data_file
            return metadata

        # Get all measurements from CSV file
        with open(data_file, 'rb') as csv_file:
            rows = csv.DictReader(csv_file)
            for row in rows:
                # Get rid
                rid = int(row['RID'])
                if rid not in metadata:
                    metadata.update({rid: {}})

                # Get scan time
                viscode = row['VISCODE']
                if viscode in metadata[rid]:
                    print log.WARNING, 'Entry already exists {0} ({1}). Skipping.'.format(rid, viscode)
                    continue
                metadata[rid].update({viscode: {}})

                # Get scan date
                scandate = datetime.datetime.strptime(row['ScanDate'], "%Y-%m-%d").date()
                metadata[rid][viscode].update({'scandate': scandate})

                # Get age
                metadata[rid][viscode].update({'AGE.scan': adni.safe_cast(row['AGE.scan'])})

                # Get factor
                metadata[rid][viscode].update({'FactorMNI': adni.safe_cast(row['FactorMNI'])})

                # Get diagnosis as numerical value
                dx = self._diagnosis_code[row['DX.scan']]
                metadata[rid][viscode].update({'DX.scan': dx})

        # Add scan time to measurements
        for rid in metadata:
            if 'bl' not in metadata[rid]:
                print log.WARNING, 'No bl scan for subject {0}!'.format(rid)
            else:
                bl_date = metadata[rid]['bl']['scandate']
                for viscode in metadata[rid]:
                    fu_date = metadata[rid][viscode]['scandate']
                    scantime = (fu_date - bl_date).days
                    metadata[rid][viscode].update({'scantime': scantime})

        # Return metadata
        print log.RESULT, 'Collected data of {0} subjects.'.format(len(metadata))
        return metadata

    ############################################################################
    #
    # update_measurements_with_biomarker_values()
    #
    ############################################################################
    def update_measurements_with_biomarker_values(self, measurements, biomarkers=None, no_regression=False):
        """ Update the metadata dictionary with the biomarker values.

        :param dict metadata: the meta data
        :param list biomarkers: the biomarkers to be read
        :param bool no_regression: do not perform age regression
        """
        biomarkers = self.get_biomarker_set() if biomarkers is None else biomarkers
        biomarker_values = self._get_biomarker_values_as_dict(measurements, biomarkers=biomarkers, no_regression=no_regression)

        # Update metadata with feature values
        for rid in measurements:
            for viscode in measurements[rid]:
                for biomarker in biomarkers:
                    try:
                        value = biomarker_values[rid][viscode][biomarker]
                        measurements[rid][viscode].update({biomarker: value})
                    except KeyError:
                        # TODO: Debug message
                        pass

        return measurements

    ############################################################################
    #
    # _get_biomarker_values_as_dict()
    #
    ############################################################################
    def _get_biomarker_values_as_dict(self, metadata, biomarkers=None, no_regression=False):
        """ Return all measurements as a dictionary.

        :param dict metadata: the meta data
        :param list biomarkers: the biomarkers to be read
        :param bool no_regression: do not perform age regression

        :return: A dictionary with the following structure:
        { <rid> : { <viscode> : { <biomarker1> : <value> }
                                { <biomarker2> : <value> } ... }
                  { <viscode> : ... }
          <rid> : ... }
        :rtype: dict
        """
        biomarkers = self.get_biomarker_set() if biomarkers is None else biomarkers

        no_regression_str = '{0} no regression'
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
                        if biomarker in adni.structure_names or biomarker in adni.structure_names_complete:
                            brain_vol = adni.safe_cast(row['Whole Brain'])
                            # factor = 1.0  # metadata[rid][viscode]['FactorMNI']
                            value *= adni.cn_mean_brain_volume / brain_vol

                        if no_regression or biomarker in adni.manifold_coordinate_names:
                            values[rid][viscode].update({biomarker: value})
                        else:
                            values[rid][viscode].update({no_regression_str.format(biomarker): value})

        if not no_regression:
            print log.INFO, 'Performing age regression...'
            for biomarker in biomarkers:
                if biomarker not in adni.manifold_coordinate_names:
                    rids = []
                    viscodes = []
                    vals = []
                    ages = []
                    for rid, visits in values.items():
                        for viscode in visits:
                            # Get age
                            try:
                                age = metadata[rid][viscode]['AGE.scan']
                            except KeyError:
                                age = None

                            if age is not None and no_regression_str.format(biomarker) in visits[viscode]:
                                value = visits[viscode][no_regression_str.format(biomarker)]
                                if value is not None:
                                    ages.append(age)
                                    vals.append(value)
                                    rids.append(rid)
                                    viscodes.append(viscode)

                    if len(ages) > 1:
                        regressed_values = self._age_regression(ages, vals)

                        for rid, viscode, value in zip(rids, viscodes, regressed_values):
                            values[rid][viscode].update({biomarker: value})

        return values

    ############################################################################
    #
    # _age_regression()
    #
    ############################################################################
    @staticmethod
    def _age_regression(timepoints, values):
        """ Perform age regression. """
        timepoints = np.array(timepoints)
        values = np.array(values)
        mean_val = np.mean(values)

        import scipy.stats as stats
        slope, intercept, _, _, _ = stats.linregress(timepoints, values)
        values = values - (timepoints * slope + intercept) + mean_val

        return values

    ############################################################################
    #
    # _select_visits()
    #
    ############################################################################
    def _select_visits(self, measurements, min_visits=0, visits=None):
        """ Select only subjects with specified visits available.

        :param dict measurements: the input measurements
        :param int min_visits: minimal number of visits that have to be available
        :param list visits: list of specified visits

        :return: the selected output measurements
        :rtype: dict
        """
        # Select subjects with minimal number of visits
        if min_visits > 0:
            print log.INFO, 'Selecting subjects with at least {0} visits...'.format(min_visits)
            valid_rids = []
            for rid in measurements:
                if len(measurements[rid]) > min_visits:
                    valid_rids.append(rid)
            measurements = {rid: value for rid, value in measurements.items() if rid in valid_rids}

        # Select subjects with specified visits
        if visits is not None:
            print log.INFO, 'Selecting subjects visits at {0}...'.format(visits)
            valid_rids = []
            for rid in measurements:
                if set(visits).issubset(set(measurements[rid].keys())):
                    valid_rids.append(rid)
            measurements = {rid: value for rid, value in measurements.items() if rid in valid_rids}

        print log.RESULT, 'Selected {0} subjects.'.format(len(measurements))
        return measurements

    ############################################################################
    #
    # _select_training_set()
    #
    ############################################################################
    def _select_training_set(self, measurements):
        """ Select only subjects converting from MCI to AD.

        :param dict measurements: the input measurements

        :return: the selected output measurements
        :rtype: dict
        """
        print log.INFO, 'Selecting converters as training set...'
        valid_rids = []
        for rid, rid_data in measurements.items():
            # Sort data according to scantime
            data_viscode = []
            data_scantime = []
            data_diagnosis = []

            for viscode, scan_data in rid_data.items():
                data_viscode.append(viscode)
                data_scantime.append(scan_data['scantime'])
                data_diagnosis.append(scan_data['DX.scan'])

            data_viscode = np.array(data_viscode)
            data_scantime = np.array(data_scantime)
            data_diagnosis = np.array(data_diagnosis)

            args = np.argsort(data_scantime)
            data_viscode = data_viscode[args]
            data_scantime = data_scantime[args]
            data_diagnosis = data_diagnosis[args]

            # Select converters with MCI as first and AD as last diagnosis
            # if data_diagnosis[-1] >= 0.25 and data_diagnosis[0] == 0.0:
            if data_diagnosis[-1] == 1.0 and 0.25 <= data_diagnosis[0] <= 0.75:
                # Mark as valid
                valid_rids.append(rid)

                # Compute time of conversion
                time_convert = None
                scantime_prev = data_scantime[0]
                for diagnosis, scantime in zip(data_diagnosis, data_scantime):
                    # if diagnosis >= 0.25:
                    if diagnosis == 1.0:
                        time_convert = scantime_prev + (scantime - scantime_prev) / 2
                        break
                    else:
                        scantime_prev = scantime

                # Update metadata with progress
                for viscode, scantime in zip(data_viscode, data_scantime):
                    measurements[rid][viscode].update({'progress': scantime - time_convert})

        # Remove non-valid rids
        metadata = {rid: value for rid, value in measurements.items() if rid in valid_rids}

        print log.RESULT, 'Selected {0} subjects.'.format(len(measurements))
        return metadata

    ############################################################################
    #
    # _select_test_set()
    #
    ############################################################################
    def _select_test_set(self, measurements):
        """ Select only non-converting AD and MCI subjects.

        :param dict measurements: the input measurements

        :return: the selected output measurements
        :rtype: dict
        """
        print log.INFO, 'Selecting non-converters as test set...'
        valid_rids = []
        for rid in measurements:
            # Sort data according to scantime
            data_scantime = []
            data_diagnosis = []

            for _, scan_data in measurements[rid].items():
                data_scantime.append(scan_data['scantime'])
                data_diagnosis.append(scan_data['DX.scan'])

            data_diagnosis = np.array(data_diagnosis)

            args = np.argsort(data_scantime)
            data_diagnosis = data_diagnosis[args]

            # Select non converters with
            if data_diagnosis[-1] == data_diagnosis[0]:
                # Mark as valid
                valid_rids.append(rid)

        # Remove non-valid rids
        metadata = {rid: value for rid, value in measurements.items() if rid in valid_rids}

        print log.RESULT, 'Selected {0} subjects.'.format(len(measurements))
        return metadata

    ############################################################################
    #
    # _select_complete_measurements()
    #
    ############################################################################
    def _select_complete_measurements(self, measurements, biomarkers=None, visits=None):
        """ Select only measurements with complete biomarker sets.

        :param dict measurements: the input measurements
        :param list biomarkers: the biomarkers that have to be available

        :return: the selected output measurements
        :rtype: dict
        """
        print log.INFO, 'Selecting subjects with complete measurements...'
        biomarkers = self.get_biomarker_set() if biomarkers is None else biomarkers

        for rid, visit in measurements.items():
            drop_rid = False
            for visit, visdata in visit.items():
                if visits is None or visit in visits:
                    for biomarker in biomarkers:
                        if biomarker not in visdata:
                            drop_rid = True
                            break
            if drop_rid:
                measurements.pop(rid)

        print log.RESULT, 'Selected {0} subjects.'.format(len(measurements))
        return measurements

    ############################################################################
    #
    # get_mmse_decline_as_dict()
    #
    ############################################################################
    @property
    def get_mmse_decline_as_dict(self):
        """ Return all a dictionary that indicates for each RID the MMSE decline
        of the subject between baseline and m24 visits.

        :return: A dictionary with the following structure:
        { <rid> : <decline>
          ... }
        :rtype: dict
        """
        measurements = self._get_metadata_as_dict(select_converters=False)
        measurements = self.update_measurements_with_biomarker_values(measurements, biomarkers=['MMSE'])

        declines = {}
        for rid in measurements:
            try:
                mmse_bl = measurements[rid]['bl']['MMSE']
                mmse_24 = measurements[rid]['m24']['MMSE']
                decline = mmse_bl - mmse_24
                declines.update({rid: decline})
            except KeyError:
                # TODO: debug message
                pass

        return declines

    ############################################################################
    #
    # get_rcd_as_dict()
    #
    ############################################################################
    @property
    def get_rcd_as_dict(self):
        """ Return all a dictionary that indicates for each RID if the subject
        is classified as RCD (rapid cognitive decline).

        :return: A dictionary with the following structure:
        { <rid> : [True|False]
          ... }
        :rtype: dict
        """
        measurements = self._get_metadata_as_dict(select_converters=False)
        measurements = self.update_measurements_with_biomarker_values(measurements, biomarkers=['MMSE'])

        rcds = {}
        for rid in measurements:
            try:
                mmse_bl = measurements[rid]['bl']['MMSE']
                mmse_24 = measurements[rid]['m24']['MMSE']
                rcd = True if (mmse_24 - mmse_bl) < -7 else False
                rcds.update({rid: rcd})
            except KeyError:
                # TODO: debug message
                pass

        return rcds


class SynthDataHandler(DataHandler):
    """
    TODO:
    """

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, args=None, iteration=None):
        """
            Initialise the right data files for the given arguments.

            :param args: command line arguments with:
            """
        super(SynthDataHandler, self).__init__()
        if args is None:
            self._biomarker_subset = None
        else:
            self._biomarker_subset = args.biomarkers_name
        self._biomarker_sets = {'synth': SynthModel.get_biomarker_names()}
        self._model_folders = {'synth': os.path.join(adni.project_folder, 'models', 'synth')}
        self._data_files = {'meta': os.path.join(adni.project_folder, 'lists/synth_data.csv'),
                            'synth': os.path.join(adni.project_folder, 'lists/synth_data.csv')}
        self._method = 'synth'
        self._volume_method = 'synth'
        if iteration is not None:
            self._iteration = iteration
        else:
            self._iteration = 0

    ############################################################################
    #
    # get_model_file()
    #
    ############################################################################
    def get_model_file(self, biomarker, iteration=None, num_samples=None, sampling=None, rate_sigma=None, run=None):
        """ Get the right model file for the given biomarker. """
        model_file = DataHandler.get_model_file(self, biomarker, iteration=iteration)
        num_samples_str = '_{0}'.format(num_samples) if num_samples is not None else ''
        sampling_str = '_{0}'.format(sampling) if sampling is not None else ''
        rate_sigma_str = '_sig{0}'.format(rate_sigma) if rate_sigma is not None and rate_sigma > 0.0 else ''
        run_str = '_{0}'.format(run) if run is not None else ''
        return model_file.replace('.csv', '{0}{1}{2}{3}.csv'.format(num_samples_str, sampling_str, rate_sigma_str, run_str))

    ############################################################################
    #
    # get_samples_file()
    #
    ############################################################################
    def get_samples_file(self, biomarker, iteration=None, num_samples=None, sampling=None, rate_sigma=None, run=None):
        """ Get the right model file for the given biomarker. """
        samples_file = DataHandler.get_samples_file(self, biomarker, iteration=iteration)
        num_samples_str = '_{0}'.format(num_samples) if num_samples is not None else ''
        sampling_str = '_{0}'.format(sampling) if sampling is not None else ''
        rate_sigma_str = '_sig{0}'.format(rate_sigma) if rate_sigma is not None and rate_sigma > 0.0 else ''
        run_str = '_{0}'.format(run) if run is not None else ''
        return samples_file.replace('.csv', '{0}{1}{2}{3}.csv'.format(num_samples_str, sampling_str, rate_sigma_str, run_str))

    ############################################################################
    #
    # _get_metadata_as_dict()
    #
    ############################################################################
    def _get_metadata_as_dict(self):
        """ Return all subjects metadata as a dictionary.

        :return: a dictionary with the following structure:
        { <rid> : { <viscode> : { DX.scan : <diagnosis> }
                                { AGE.scan : <age in years> }
                                { scandate : <date of scan> }
                                { scantime : <days after bl> }
                                { progress : <days relative to conversion> } }
                  { <viscode> : ... }
          <rid> : ... }
        :rtype: dict
        """
        print log.INFO, 'Collecting metadata...'
        metadata = {}

        data_file = self.get_data_file(method='meta')
        if not os.path.isfile(data_file):
            print log.ERROR, 'Data file dies not exist:', data_file
            return metadata

        #
        # Get all measurements from CSV file
        with open(data_file, 'rb') as csv_file:
            rows = csv.DictReader(csv_file)
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
                dx = self._diagnosis_code[row['DX.scan']]
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
        """ Return all measurements as a dictionary.

        :param dict metadata: the meta data
        :param list biomarkers: the biomarkers to be read
        :param bool no_regression: do not perform age regression

        :return: a dictionary with the following structure:
        { <rid> : { <viscode> : { <biomarker1> : <value> }
                                { <biomarker2> : <value> } ... }
                  { <viscode> : ... }
          <rid> : ... }
        :rtype: dict
        """
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
    # _select_training_set()
    #
    ############################################################################
    def _select_training_set(self, measurements):
        """ Select only subjects converting from MCI to AD.

        :param dict measurements: the input measurements

        :return: the selected output measurements
        :rtype: dict
        """
        return measurements

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
