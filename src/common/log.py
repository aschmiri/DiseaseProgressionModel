#! /usr/bin/env python2.7
from sys import stdout

################################################################################
#
# Colour output
#
################################################################################
SKIP = '\033[95mSKIP:\033[0m' if stdout.isatty() else 'SKIP:'
INFO = '\033[94mINFO:\033[0m' if stdout.isatty() else 'INFO:'
RESULT = '\033[92mRESULT:\033[0m' if stdout.isatty() else 'RESULT:'
WARNING = '\033[93mWARNING:\033[0m' if stdout.isatty() else 'WARNING:'
ERROR = '\033[91mERROR:\033[0m' if stdout.isatty() else 'ERROR:'
