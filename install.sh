#!/bin/bash
virtualenv env
. env/bin/activate
pip install numpy
pip install cython
pip install numexpr
pip install tables
