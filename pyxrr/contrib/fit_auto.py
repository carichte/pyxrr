#!/usr/bin/env python

# Copyright 2012 by Carsten Richter
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
import sys

if os.path.isfile("../pyxrr/xrr.so"):
    sys.path.insert(0, os.path.abspath("../"))

import pickle
import pyxrr
import pylab

"""
    This is a rather advanced example for using pyxrr to fit measured data automaticely.
    
    
    see pyxrr.multilayer.__init__.__doc__ for explanation.
"""

FITTYPE = "log"
VERBOSE = 2
PENALTY = 1.
DATABASE_f1f2 = "Henke" # can be BrennanCowan, Chantler, CromerLiberman, EPDL97, Henke, Sasaki, Windt
                        # For more details see Databases on http://www.esrf.eu/computing/scientific/dabax
                        # and http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/







################################################################################


ROOT="WSi_stack2"
SAMPLENAME=os.path.join("samples", ROOT + ".param")

sample=pyxrr.multilayer(SAMPLENAME, DB_Table = DATABASE_f1f2, verbose=VERBOSE, penalty=PENALTY, fittype=FITTYPE)



### SET COUPLED AND VARIABLE PARAMETERS
fit_variables=[]

for num in range(4,33,2): sample.coupled_vars["rho_%i"%num] = "rho_2" # Alle Wolfram-Schichten haben gleiche Dichte
fit_variables.append("rho_2")

for num in range(5,34,2): sample.coupled_vars["rho_%i"%num] = "rho_3" # Alle Silizium-Spacer Schichten haben gleiche Dichte
fit_variables.append("rho_3")

for num in range(14): 
    sample.coupled_vars["sigma_%i"%(33-num)] = "sigma_%i"%(15-num) # Ein paar Grenzflaechen sind auch noch equivalent (Da 2 Gruppen)
    fit_variables.append("sigma_%i"%(15-num))
    sample.coupled_vars["d_%i"%(33-num)] = "d_%i"%(15-num) # Gleiches gilt fuer die Schichtdicken
    fit_variables.append("d_%i"%(15-num))
    

#print(sample.coupled_vars) # Kontrolle
#print(fit_variables)

#sample.couple()
#print(sample.print_parameter())

# noch ein paar Parameter freischalten
fit_variables+=["sigma_1", "sigma_0", "d_1", "rho_1"]
fit_variables+=["sigma_19", "sigma_18", "sigma_17", "sigma_16", "d_19", "d_18", "d_17", "d_16"]
fit_variables+=["offset0", "offset1", "offset2"]
fit_variables+=["background0", "background1", "background2"]




### START LEAST-SQUARES FIT:
sample.fit("leastsq", fit_variables)



### SAVE RESULTS
fn=open(os.path.join("results", ROOT + ".save"), "w")
pickle.dump(sample.parameters, fn)
fn.close()
print "Saved Parameters to " + os.path.join("results", ROOT + ".save")

fn=open(os.path.join("results", ROOT + ".txt"), "w")
fn.write(sample.print_parameter())
try: fn.write("\nError: " + str((sample.err**2).sum()))
except: fn.write("\nError: " + str((sample.residuals({})**2).sum()))
fn.close()
print "Saved Parameters (ASCII) to " + os.path.join("results", ROOT + ".txt")

for i_M in range(sample.number_of_measurements): 
    pylab.savetxt(os.path.join("results", ROOT + "_M%i.dat"%i_M), pylab.vstack((sample.measured_data[i_M][:,0], sample.reflectogram(i_M=i_M), sample.measured_data[i_M][:,1])).T)

