#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

from glasson import _VERSION


setup(
    name='glasson',
    version=_VERSION,
    description='glasson is a binary and compressed file format for storing Hi-C contact maps',
    author='Matthias Blum',
    py_module='glasson',
    scripts=['bin/glatools'],
    url='https://github.com/matthiasblum/glasson',
    license='Public Domain'
)
