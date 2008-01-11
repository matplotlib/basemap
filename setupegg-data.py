"""
Poor man's setuptools script...
"""

from setuptools import setup
execfile('setup-data.py',
         {'additional_params' :
         {'namespace_packages' : ['mpl_toolkits']}})
