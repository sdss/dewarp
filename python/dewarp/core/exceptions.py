# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date:   2017-12-05 12:01:21
# @Last modified by:   Brian Cherinka
# @Last Modified time: 2017-12-05 12:19:32

from __future__ import print_function, division, absolute_import


class DewarpError(Exception):
    """A custom core Dewarp exception"""

    def __init__(self, message=None):

        message = 'There has been an error' \
            if not message else message

        super(DewarpError, self).__init__(message)


class DewarpNotImplemented(DewarpError):
    """A custom exception for not yet implemented features."""

    def __init__(self, message=None):

        message = 'This feature is not implemented yet.' \
            if not message else message

        super(DewarpNotImplemented, self).__init__(message)


class DewarpAPIError(DewarpError):
    """A custom exception for API errors"""

    def __init__(self, message=None):
        if not message:
            message = 'Error with Http Response from Dewarp API'
        else:
            message = 'Http response error from Dewarp API. {0}'.format(message)

        super(DewarpAPIError, self).__init__(message)


class DewarpApiAuthError(DewarpAPIError):
    """A custom exception for API authentication errors"""
    pass


class DewarpMissingDependency(DewarpError):
    """A custom exception for missing dependencies."""
    pass


class DewarpWarning(Warning):
    """Base warning for Dewarp."""


class DewarpUserWarning(UserWarning, DewarpWarning):
    """The primary warning class."""
    pass


class DewarpSkippedTestWarning(DewarpUserWarning):
    """A warning for when a test is skipped."""
    pass


class DewarpDeprecationWarning(DewarpUserWarning):
    """A warning for deprecated features."""
    pass
