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


"""
    EXAMPLE for using pyxrr
    this batchfile searches for optimal gamma parameters for sample def06
    gamma = d_absorber/d_spacer
"""
import os.path
import sys

if os.path.isfile("../pyxrr/xrr.so"):
    sys.path.insert(0, os.path.abspath("../"))


from pyxrr import *
from scipy.optimize import fmin

sample=multilayer("samples/def06.param")


d=sample.parameters["d_2"]+sample.parameters["d_3"]


results=np.zeros(4)


def Rmax_inv(v):
    [theta, gamma]=v
    sample.parameters["d_2"]=gamma*d
    sample.parameters["d_3"]=(1.-gamma)*d
    thisRmax=sample.reflectogram(np.array([theta]))
    return 1./thisRmax
    

for energy in np.arange(1,24.01,0.1):
    th_max=np.arcsin(12.398/energy/(2*d))*180/np.pi
    sample.parameters["energy0"]=energy
    ([th_opt, gamma_opt], Rmax, iters, funcalls, warnflag) = fmin(Rmax_inv, [th_max, 0.3], xtol=1e-8, full_output=1)
    Rmax=1/Rmax
    print energy, gamma_opt, Rmax
    results=np.vstack((results,np.array([energy,th_opt,gamma_opt,Rmax])))

np.savetxt("gammamax_def06.txt",results)

