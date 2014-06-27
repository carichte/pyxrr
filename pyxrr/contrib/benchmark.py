#!/usr/bin/env python

# Copyright 2013 by Carsten Richter
# Contact: carsten.richter@desy.de and
#          carsten.richter@physik.tu-freiberg.de
#
# This file is part of pyxrr.
#
# pyxrr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyxrr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyxrr inside the 'copying.txt' file.
# If not, see <http://www.gnu.org/licenses/>.


import os
#os.environ["OPENBLAS_MAIN_FREE"] = '1'
import sys
import warnings
import pickle


sys.path.insert(0, os.path.abspath(os.pardir)) # Path to pyxrr if it has not been installed as package
try:
    import pyxrr
except:
    sys.path.pop(0) # not there
    import pyxrr


warnings.filterwarnings("ignore")


"""
    This is a simple batch file for using pyxrr to fit measured data 
    interactively. It can be seen as an example
    
    see pyxrr.multilayer.__init__.__doc__ for an explanation.
"""

FITTYPE = "log"
VERBOSE = 2
PENALTY = 1.
DATABASE_f1f2 = "Henke" # can be BrennanCowan, Chantler, CromerLiberman, EPDL97, Henke, Sasaki, Windt
                        # For more details see Databases on http://www.esrf.eu/computing/scientific/dabax
                        # and http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/






################################################################################




if not os.path.isdir("tmp"):
    os.makedirs("tmp")

SAMPLENAME = "test24.param"
ROOT = SAMPLENAME.rpartition(".")[0]
if not ROOT:
    ROOT = SAMPLENAME


sample=pyxrr.multilayer(os.path.join("samples",SAMPLENAME), 
                        DB_Table = DATABASE_f1f2, verbose=VERBOSE, 
                        penalty=PENALTY, fittype=FITTYPE)


N=1000
t0 = pyxrr.time.time()

for i in xrange(N):
    R = sample.reflectogram()

dt = pyxrr.time.time() - t0
print("Time per function call: %g ms"%(1000.*dt/N))
