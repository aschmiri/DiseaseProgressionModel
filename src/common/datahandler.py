"""
A class to provide reading and writing data for progression modelling.

:author:     Alexander Schmidt-Richberg
:copyright:  2014 Imperial College London. All rights reserved.
:contact:    a.schmidt-richberg@imperial.ac.uk
"""
import os.path
import datetime
import csv
import ConfigParser
import numpy as np
import scipy.stats as stats
from common import log as log
from common.synthmodel import SynthModel


############################################################################
#
# DataHandler()
#
############################################################################
class DataHandler(object):
    """
    TODO: classdocs
    """
    class Config(object):
        class ListConfigParser(ConfigParser.SafeConfigParser):
            def getlist(self, section, option):
                list_string = self.get(section, option)
                return [e.strip() for e in list_string.split(',')]

            def getint(self, section, option):
                value = self.get(section, option)
                try:
                    return int(value)
                except:
                    return None

            def getfloat(self, section, option):
                value = self.get(section, option)
                try:
                    return float(value)
                except:
                    return None

        def __init__(self):
            # Load config file
            config_path = os.path.dirname(os.path.relpath(__file__))
            config_file = os.path.join(config_path, '..', 'configure.ini')
            config = self.ListConfigParser()
            config.read(config_file)

            # Load default data and project folders
            self.project_folder = config.get('DEFAULT', 'project_folder')
            self.data_folder = config.get('DEFAULT', 'data_folder')
            self.models_folder = config.get('DEFAULT', 'models_folder')
            self.eval_folder = config.get('DEFAULT', 'eval_folder')

            # Load biomarker configurations
            self.biomarker_sets = config.getlist('biomarkers', 'biomarker_sets')
            self.combined_biomarker_sets = config.getlist('biomarkers', 'combined_sets')
            self.biomarker_names = {'synth': SynthModel.get_biomarker_names()}
            self.model_folders = {'synth': os.path.join(self.models_folder, 'synth')}
            self.data_files = {'synth': os.path.join(self.data_folder, 'synth_data.csv'),
                               'meta': config.get('biomarkers', 'meta_data_file')}
            self.age_regression = {'synth': False}
            self.biomarker_units = {'synth': 'Score'}

            for biomarker_set in self.biomarker_sets:
                data_file = config.get(biomarker_set, 'data_file')
                self.data_files.update({biomarker_set: data_file})

                biomarker_names = config.getlist(biomarker_set, 'biomarker_names')
                self.biomarker_names.update({biomarker_set: biomarker_names})

                age_regression = config.get(biomarker_set, 'age_regression')
                age_regression = True if age_regression in ['Yes', 'yes', 'True', 'true'] else False
                self.age_regression.update({biomarker_set: age_regression})

                biomarker_unit = config.get(biomarker_set, 'biomarker_unit')
                self.biomarker_units.update({biomarker_set: biomarker_unit})

            for combined_set in self.combined_biomarker_sets:
                biomarker_names = []
                if config.has_option(combined_set, 'biomarker_sets'):
                    biomarker_sets = config.getlist(combined_set, 'biomarker_sets')
                    for biomarker_set in biomarker_sets:
                        biomarker_names += config.getlist(biomarker_set, 'biomarker_names')
                    self.biomarker_names.update({combined_set: biomarker_names})
                if config.has_option(combined_set, 'biomarker_names'):
                    biomarker_names += config.getlist(combined_set, 'biomarker_names')
                if len(biomarker_names) == 0:
                    print log.ERROR, 'Failed to read set {0}'.format(combined_set)
                else:
                    self.biomarker_names.update({combined_set: biomarker_names})

            # Load VGAM configurations
            self.model_offset = config.getint('VGAM', 'model_offset')
            self.vgam_degrees_of_freedom = {'default': config.getint('VGAM', 'degrees_of_freedom')}
            self.vgam_zero = {'default': config.getint('VGAM', 'zero')}
            for biomarker_set in self.biomarker_sets:
                # Load specific settings for biomarker set
                if config.has_option(biomarker_set, 'degrees_of_freedom'):
                    dof = config.getint(biomarker_set, 'degrees_of_freedom')
                    self.vgam_degrees_of_freedom.update({biomarker_set: dof})
                if config.has_option(biomarker_set, 'zero'):
                    zero = config.getint(biomarker_set, 'zero')
                    self.vgam_zero.update({biomarker_set: zero})

                # Load specific settings for biomarker
                for biomarker in self.biomarker_names[biomarker_set]:
                    if config.has_section(biomarker):
                        if config.has_option(biomarker, 'degrees_of_freedom'):
                            dof = config.getint(biomarker, 'degrees_of_freedom')
                            self.vgam_degrees_of_freedom .update({biomarker: dof})
                        if config.has_option(biomarker, 'zero'):
                            zero = config.getint(biomarker, 'zero')
                            self.vgam_zero.update({biomarker: zero})

    _conf = Config()
    _method = None
    _biomarker_subset = None

    ############################################################################
    #
    # get_data_handler()
    #
    ############################################################################
    @staticmethod
    def get_data_handler(method=None, biomarkers=None, phase=None):
        if method == 'synth':
            return SynthDataHandler(biomarkers=biomarkers)
        else:
            return ClinicalDataHandler(method=method, biomarkers=biomarkers, phase=phase)

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, biomarkers=None):
        """
        Initialise the right data files for the given arguments.

        :param args: command line arguments with:
        """
        self._biomarker_subset = biomarkers
        self._diagnosis_code = {'CN': 0.0, 'EMCI': 0.25, 'MCI': 0.5, 'LMCI': 0.75, 'AD': 1.0}

    ############################################################################
    #
    # get_method_choices()
    #
    ############################################################################
    @classmethod
    def get_method_choices(cls):
        return cls._conf.biomarker_names.keys()

    ############################################################################
    #
    # get_method_choices()
    #
    ############################################################################
    @classmethod
    def get_phase_choices(cls):
        return ['cnmci', 'mciad', 'joint']

    ############################################################################
    #
    # safe_cast
    #
    ############################################################################
    @staticmethod
    def safe_cast(in_value, cast_type=float):
        try:
            return cast_type(in_value)
        except ValueError:
            return None
        except TypeError:
            return None

    ############################################################################
    #
    # make_dir
    #
    ############################################################################
    @staticmethod
    def make_dir(*dirs):
        directory = os.path.join(*dirs)
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory

    ############################################################################
    #
    # get_project_folder()
    #
    ############################################################################
    def get_project_folder(self):
        """ Get the project folder."""
        return self.make_dir(self._conf.project_folder)

    ############################################################################
    #
    # get_data_folder()
    #
    ############################################################################
    def get_data_folder(self):
        """ Get the data folder."""
        return self.make_dir(self._conf.data_folder)

    ############################################################################
    #
    # _get_model_folder_for()
    #
    ############################################################################
    def get_model_folder(self, biomarker):
        """ Get the right method folder for the given biomarker. """
        return self.make_dir(self._conf.models_folder, self._get_method_for_biomarker(biomarker))

    ############################################################################
    #
    # get_eval_folder()
    #
    ############################################################################
    def get_eval_folder(self):
        """ Get the evaluation folder."""
        return self.make_dir(self._conf.eval_folder)

    ############################################################################
    #
    # get_model_file()
    #
    ############################################################################
    def get_model_file(self, biomarker):
        """ Get the right model file for the given biomarker. """
        model_folder = self.get_model_folder(biomarker)
        return os.path.join(model_folder, biomarker.replace(' ', '_') + '_model.csv')

    ############################################################################
    #
    # get_samples_file()
    #
    ############################################################################
    def get_samples_file(self, biomarker):
        """ Get the right model file for the given biomarker. """
        model_folder = self.get_model_folder(biomarker)
        return os.path.join(model_folder, biomarker.replace(' ', '_') + '_samples.csv')

    ############################################################################
    #
    # _get_data_file_for_biomarker()
    #
    ############################################################################
    def _get_data_file_for_biomarker(self, biomarker):
        """ Get the right data file for the given biomarker. """
        method = self._get_method_for_biomarker(biomarker)
        return self._get_data_file_for_method(method)

    ############################################################################
    #
    # _get_data_file_for_method()
    #
    ############################################################################
    def _get_data_file_for_method(self, method):
        """ Get the right data file for the given method. """
        return self._conf.data_files[method]

    ############################################################################
    #
    # _get_method_for_biomarker()
    #
    ############################################################################
    @classmethod
    def _get_method_for_biomarker(cls, biomarker):
        """ Get the right method for the given biomarker."""
        # Test names of all biomarker sets
        for biomarker_set in cls._conf.biomarker_sets:
            if biomarker in cls._conf.biomarker_names[biomarker_set]:
                return biomarker_set

        # Test names of synth models
        if biomarker in cls._conf.biomarker_names['synth']:
                return 'synth'

        # Return error if not found
        print log.ERROR, 'Method cannot be determined for {0}!'.format(biomarker)
        return None

    ############################################################################
    #
    # get_biomarker_unit()
    #
    ############################################################################
    @classmethod
    def get_biomarker_unit(cls, biomarker):
        """ Get the right unit for the given biomarker, mainly required for plotting. """
        return cls._conf.biomarker_units[cls._get_method_for_biomarker(biomarker)]

    ############################################################################
    #
    # get_all_biomarker_names()
    #
    ############################################################################
    @classmethod
    def get_all_biomarker_names(cls):
        """ Get the all biomarker names."""
        biomarker_names = []
        for biomarker_set in cls._conf.biomarker_sets:
            biomarker_names += cls._conf.biomarker_names[biomarker_set]
        return biomarker_names

    ############################################################################
    #
    # get_biomarker_names()
    #
    ############################################################################
    def get_biomarker_names(self, method=None):
        """ Get the right set of biomarkers for the current configuration."""
        if method is not None:
            return self._conf.biomarker_names[method]
        elif self._biomarker_subset is not None:
            return self._biomarker_subset
        else:
            return self._conf.biomarker_names[self._method]

    ############################################################################
    #
    # get_vgam_degrees_of_freedom()
    #
    ############################################################################
    def get_vgam_degrees_of_freedom(self, biomarker=None):
        """ Get the degrees of freedom for the given biomarker.
        If conflicting settings are given, the following order is considered:
        biomarker-specific > set-specific > default definition
        """
        if biomarker is None:
            return self._conf.vgam_degrees_of_freedom['default']
        elif biomarker in self._conf.vgam_degrees_of_freedom:
            return self._conf.vgam_degrees_of_freedom[biomarker]
        elif self._get_method_for_biomarker(biomarker) in self._conf.vgam_degrees_of_freedom:
            return self._conf.vgam_degrees_of_freedom[self._get_method_for_biomarker(biomarker)]
        else:
            return self._conf.vgam_degrees_of_freedom['default']

    ############################################################################
    #
    # get_vgam_zero()
    #
    ############################################################################
    def get_vgam_zero(self, biomarker=None):
        """ Get the zero for the given biomarker.
        If conflicting settings are given, the following order is considered:
        biomarker-specific > set-specific > default definition
        """
        if biomarker is None:
            return self._conf.vgam_zero['default']
        elif biomarker in self._conf.vgam_zero:
            return self._conf.vgam_zero[biomarker]
        elif self._get_method_for_biomarker(biomarker) in self._conf.vgam_zero:
            return self._conf.vgam_zero[self._get_method_for_biomarker(biomarker)]
        else:
            return self._conf.vgam_zero['default']

    ############################################################################
    #
    # get_measurements_as_dict()
    #
    ############################################################################
    def get_measurements_as_dict(self, **kwargs):
        """ Return all subjects measurements as a dictionary.

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
        return {}


############################################################################
#
# ClinicalDataHandler()
#
############################################################################
class ClinicalDataHandler(DataHandler):
    """
    TODO: classdocs
    """

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, method=None, biomarkers=None, phase=None):
        """
        Initialise the right data files for the given arguments.

        :param str args.method: one of the methods defined in the config file
        :param list args.biomarkers: subset of biomarkers defined in the config file
        """
        super(ClinicalDataHandler, self).__init__(biomarkers)

        # Set method
        if method is None:
            self._method = 'all'
        else:
            self._method = method

        if phase is None:
            self._phase = 'mciad'
        else:
            self._phase = phase

    ############################################################################
    #
    # get_model_folder()
    #
    ############################################################################
    def get_model_folder(self, biomarker):
        """ Get the right method folder for the given biomarker. """
        return self.make_dir(self._conf.models_folder,
                             self._phase,
                             self._get_method_for_biomarker(biomarker))

    ############################################################################
    #
    # get_eval_folder()
    #
    ############################################################################
    def get_eval_folder(self):
        """ Get the evaluation folder."""
        return self.make_dir(self._conf.eval_folder,
                             self._phase)

    ############################################################################
    #
    # get_model_offset()
    #
    ############################################################################
    def get_model_offset(self):
        return self._conf.model_offset

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
        measurements = self._update_measurements_with_biomarker_values(measurements,
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

        data_file = self._get_data_file_for_method('meta')
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
                metadata[rid][viscode].update({'AGE.scan': self.safe_cast(row['AGE.scan'])})

                # Get factor
                # metadata[rid][viscode].update({'FactorMNI': self.safe_cast(row['FactorMNI'])})

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
    # _update_measurements_with_biomarker_values()
    #
    ############################################################################
    def _update_measurements_with_biomarker_values(self, measurements, biomarkers=None, no_regression=False):
        """ Update the metadata dictionary with the biomarker values.

        :param dict metadata: the meta data
        :param list biomarkers: the biomarkers to be read
        :param bool no_regression: do not perform age regression
        """
        biomarkers = self.get_biomarker_names() if biomarkers is None else biomarkers
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
        biomarkers = self.get_biomarker_names() if biomarkers is None else biomarkers

        no_regression_str = '{0} no regression'
        print log.INFO, 'Collecting biomarker values...'
        values = {}

        for biomarker in biomarkers:
            data_file = self._get_data_file_for_biomarker(biomarker)
            if data_file is None:
                print log.ERROR, 'Data file not found:', data_file
                continue

            age_regression = self._conf.age_regression[self._get_method_for_biomarker(biomarker)]

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
                        value = self.safe_cast(row[biomarker])
                        if value is not None:
                            if age_regression and not no_regression:
                                values[rid][viscode].update({no_regression_str.format(biomarker): value})
                            else:
                                values[rid][viscode].update({biomarker: value})

        for biomarker in biomarkers:
            if not no_regression and self._conf.age_regression[self._get_method_for_biomarker(biomarker)]:
                print log.INFO, 'Performing age regression for {0}...'.format(biomarker)
                rids = []
                viscodes = []
                vals = []
                ages = []
                diagnoses = []
                for rid, visits in values.items():
                    for viscode in visits:
                        # Get age
                        try:
                            age = metadata[rid][viscode]['AGE.scan']
                            diagnosis = metadata[rid][viscode]['DX.scan']
                        except KeyError:
                            age = None

                        if age is not None and no_regression_str.format(biomarker) in visits[viscode]:
                            value = visits[viscode][no_regression_str.format(biomarker)]
                            if value is not None:
                                ages.append(age)
                                vals.append(value)
                                rids.append(rid)
                                viscodes.append(viscode)
                                diagnoses.append(diagnosis)

                if len(ages) > 1:
                    regressed_values = self._age_regression(ages, vals, viscodes, diagnoses)

                    for rid, viscode, value in zip(rids, viscodes, regressed_values):
                        values[rid][viscode].update({biomarker: value})

        return values

    ############################################################################
    #
    # _age_regression()
    #
    ############################################################################
    def _age_regression(self, timepoints, values, viscodes, diagnoses):
        """ Perform age regression. """
        # Cast input lists to numpy arrays
        timepoints = np.array(timepoints)
        values = np.array(values)
        diagnoses = np.array(diagnoses)
        viscodes = np.array(viscodes)

        # Compute linear regression model based on baseline values of CN subjects
        cn_indices = np.where((diagnoses == self._diagnosis_code['CN']) & (viscodes == 'bl'))
        cn_timepoints = timepoints[cn_indices]
        cn_values = values[cn_indices]
        cn_mean = np.mean(cn_values)
        slope, intercept, _, _, _ = stats.linregress(cn_timepoints, cn_values)

        # Transform values with progression line
        values = values - (timepoints * slope + intercept) + cn_mean
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

            converter_cn_mci = self._diagnosis_is_cn(data_diagnosis[0]) and self._diagnosis_is_mci(data_diagnosis[-1])
            converter_cn_ad = self._diagnosis_is_cn(data_diagnosis[0]) and self._diagnosis_is_ad(data_diagnosis[-1])
            converter_mci_ad = self._diagnosis_is_mci(data_diagnosis[0]) and self._diagnosis_is_ad(data_diagnosis[-1])
            converter = converter_cn_mci or converter_mci_ad

            # Select converters with MCI as first and AD as last diagnosis
            if ((self._phase == 'cnmci' and (converter_cn_mci or converter_cn_ad)) or
                    (self._phase == 'mciad' and converter_mci_ad) or
                    (self._phase == 'joint' and converter)):
                # Mark as valid
                valid_rids.append(rid)

                # Compute time of conversion
                time_convert = None
                scantime_prev = data_scantime[0]
                for diagnosis, scantime in zip(data_diagnosis, data_scantime):
                    if ((self._phase == 'cnmci' and
                         (self._diagnosis_is_mci(diagnosis) or self._diagnosis_is_ad(diagnosis))) or
                        (self._phase == 'mciad' and
                         self._diagnosis_is_ad(diagnosis))):
                        time_convert = 0.5 * (scantime + scantime_prev)
                        break
                    elif (self._phase == 'joint' and converter_cn_mci and
                          self._diagnosis_is_mci(diagnosis)):
                        time_convert = 0.5 * (scantime + scantime_prev)
                        break
                    elif (self._phase == 'joint' and converter_mci_ad and
                          self._diagnosis_is_ad(diagnosis)):
                        time_convert = 0.5 * (scantime + scantime_prev) - self._conf.model_offset
                        break
                    else:
                        scantime_prev = scantime

                # Update metadata with progress
                for viscode, scantime in zip(data_viscode, data_scantime):
                    measurements[rid][viscode].update({'progress': scantime - time_convert})

        # Remove non-valid rids
        measurements = {rid: value for rid, value in measurements.items() if rid in valid_rids}

        print log.RESULT, 'Selected {0} subjects.'.format(len(measurements))
        return measurements

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
        measurements = {rid: value for rid, value in measurements.items() if rid in valid_rids}

        print log.RESULT, 'Selected {0} subjects.'.format(len(measurements))
        return measurements

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
        biomarkers = self.get_biomarker_names() if biomarkers is None else biomarkers

        for rid, visit in measurements.items():
            drop_rid = False
            for viscode, visdata in visit.items():
                if visits is None or viscode in visits:
                    if not set(biomarkers).issubset(set(visdata)):
                        drop_rid = True
            if drop_rid:
                measurements.pop(rid)

        print log.RESULT, 'Selected {0} subjects.'.format(len(measurements))
        return measurements

    ############################################################################
    #
    # _diagnosis_is_cn()
    #
    ############################################################################
    def _diagnosis_is_cn(self, diagnosis):
        return self._diagnosis_code['CN'] == diagnosis

    ############################################################################
    #
    # _diagnosis_is_mci()
    #
    ############################################################################
    def _diagnosis_is_mci(self, diagnosis):
        return self._diagnosis_code['EMCI'] <= diagnosis <= self._diagnosis_code['LMCI']

    ############################################################################
    #
    # _diagnosis_is_ad()
    #
    ############################################################################
    def _diagnosis_is_ad(self, diagnosis):
        return self._diagnosis_code['AD'] == diagnosis


class SynthDataHandler(DataHandler):
    """
    TODO:
    """

    ############################################################################
    #
    # __init__()
    #
    ############################################################################
    def __init__(self, biomarkers=None):
        """
        Initialise the right data files for the given arguments.

        :param args: command line arguments with:
        """
        super(SynthDataHandler, self).__init__(biomarkers)
        self._method = 'synth'

    ############################################################################
    #
    # get_model_file()
    #
    ############################################################################
    def get_model_file(self, biomarker, num_samples=None, sampling=None,
                       rate_sigma=None, conversion_sigma=None, run=None):
        """ Get the right model file for the given biomarker. """
        model_file = DataHandler.get_model_file(biomarker)
        num_samples_str = '_{0}'.format(num_samples) if num_samples is not None else ''
        sampling_str = '_{0}'.format(sampling) if sampling is not None else ''
        rate_sigma_str = '_sig{0}'.format(rate_sigma) if rate_sigma is not None and rate_sigma > 0.0 else ''
        conv_sigma_str = '_csig{0}'.format(conversion_sigma) if conversion_sigma is not None and conversion_sigma > 0.0 else ''
        run_str = '_{0}'.format(run) if run is not None else ''
        return model_file.replace('.csv', '{0}{1}{2}{3}{4}.csv'.format(
            num_samples_str, sampling_str, rate_sigma_str, conv_sigma_str, run_str))

    ############################################################################
    #
    # get_samples_file()
    #
    ############################################################################
    def get_samples_file(self, biomarker, num_samples=None, sampling=None,
                         rate_sigma=None, conversion_sigma=None, run=None):
        """ Get the right model file for the given biomarker. """
        samples_file = DataHandler.get_samples_file(biomarker)
        num_samples_str = '_{0}'.format(num_samples) if num_samples is not None else ''
        sampling_str = '_{0}'.format(sampling) if sampling is not None else ''
        rate_sigma_str = '_sig{0}'.format(rate_sigma) if rate_sigma is not None and rate_sigma > 0.0 else ''
        conv_sigma_str = '_csig{0}'.format(conversion_sigma) if conversion_sigma is not None and conversion_sigma > 0.0 else ''
        run_str = '_{0}'.format(run) if run is not None else ''
        return samples_file.replace('.csv', '{0}{1}{2}{3}{4}.csv'.format(
            num_samples_str, sampling_str, rate_sigma_str, conv_sigma_str, run_str))

    ############################################################################
    #
    # get_measurements_as_dict()
    #
    ############################################################################
    def get_measurements_as_dict(self, biomarkers=None, **kwargs):
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
        biomarkers = self.get_biomarker_names() if biomarkers is None else biomarkers
        measurements = {}

        data_file = self._conf.data_files['synth']
        if not os.path.isfile(data_file):
            print log.ERROR, 'Data file dies not exist:', data_file
            return measurements

        #
        # Get all measurements from CSV file
        with open(data_file, 'rb') as csv_file:
            rows = csv.DictReader(csv_file)
            for row in rows:
                # Get rid
                rid = int(row['RID'])
                if rid not in measurements:
                    measurements.update({rid: {}})

                # Get scan time
                viscode = row['VISCODE']
                if viscode in measurements[rid]:
                    print log.WARNING, 'Entry already exists {0} ({1}). Skipping.'.format(rid, viscode)
                    break
                measurements[rid].update({viscode: {}})

                # Get diagnosis as numerical value
                dx = self._diagnosis_code[row['DX.scan']]
                measurements[rid][viscode].update({'DX.scan': dx})

                # Get progress
                progress = row['progress']
                measurements[rid][viscode].update({'progress': progress})

                for biomarker in biomarkers:
                    value = DataHandler.safe_cast(row[biomarker])
                    measurements[rid][viscode].update({biomarker: value})

        return measurements
