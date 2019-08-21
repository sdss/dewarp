# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      August 21, 2019
# @Filename:  __init__.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from .utils import get_config
from .utils import get_logger

NAME = 'dewarp'

# Loads config
config = get_config(NAME)

# Inits the logging system. Only shell logging, and exception and warning catching.
# File logging can be started by calling log.start_file_logger(path).
log = get_logger(NAME)

__version__ = '0.0.0'
