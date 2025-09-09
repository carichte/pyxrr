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
import time
import json
import numpy as np
import lmfit
import collections
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d

from .xrr import reflectivity as _reflectivity
from . import xray_interactions as xi
from . import measurement, structure


try:
    import pandas as pd
    PANDAS = True
except:
    PANDAS = False

FFT = np.fft.rfft


class pyxrrError(Exception):
    def __init__(self, value, errmsg="", identifier=None):
        self.value = value
        self.errmsg = errmsg
    def __str__(self):
        return self.value + lsep + "Message: %s" %self.errmsg




class Model(object):
    """
        Main XRR calc/fit class
    """
    def __init__(self,
                 stack,
                 measurements,
                 table="Henke", 
                 numthreads=0,
                 verbose=True,
                 benchmark=False):
        """
            multilayer object as defined in the given 'SampleFile'.
            
            Optional inputs:
            
             - table (string) -  set of data to use for dispersion correction
                                    coefficients. Can be:
                - 'BrennanCowan' (see http://skuld.bmsc.washington.edu/scatter/)
                - 'Chantler' (see http://www.nist.gov/pml/data/ffast/index.cfm)
                - 'CromerLiberman'
                - 'EPDL97'
                - 'Henke' (see http://www-cxro.lbl.gov/optical_constants/)
                - 'Sasaki'
                - 'Windt' (also used by the Program IMD, see f.e. 
                            http://www.rxollc.com/idl/)
                
                For more details see Databases on 
                http://www.esrf.eu/computing/scientific/dabax and 
                http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/
        """
        if isinstance(measurements, measurement.Measurement):
            measurements = [measurements]

        self.stack = stack
        stack.update()
        self.measurements = dict((m.id, m) for m in measurements)
        self.table = table
        self.numthreads = numthreads
        self.update()

    def update(self):
        """
            Creates a new `Parameters` object that contains the
            parameters describing both the stack and the measurement.
        """
        params = self.stack.params
        for m in self.measurements.values():
            params.add_many(*m.get_params())
        self.params = params

    def update_params(self, new_params, values_only=True):
        if isinstance(new_params, lmfit.Parameters):
            for key in new_params:
                self.params[key].value = new_params[key].value
                if not values_only:
                    self.params[key].min = new_params[key].min
                    self.params[key].max = new_params[key].max
                    self.params[key].vary = new_params[key].vary
                    self.params[key].expr = new_params[key].expr
        else:
            for key in new_params:
                self.params[key].value = new_params[key]


    def dump(self, path, overwrite=False):
        if os.path.isfile(path) and not overwrite:
            raise IOError("File already exists. Use `overwrite` argument.")
        with open(path, 'w') as f:
            self.params.dump(f)

    def load(self, path):
        params = lmfit.Parameters()
        with open(path, 'r') as f:
            params.load(f)

        self.update_params(params, False)

    def fetch_optical_constants(self, energy, table=None):
        if table is None:
            table = self.table
        stack = self.stack
        return xi.get_optical_constants(stack._compositions, energy, stack._densities)



    def reflectivity(self, x=None, idm=0, **new_params):
        """
            Calculates a reflectivity curve for a given theta range.
            The theta array has to be equally spaced data and sorted.
            `idm` specifies the measurement from which to take parameters like 
            energy, offset, polarization etc.
        """
        if new_params:
            self.update_params(new_params)

        m = self.measurements[idm]
        stack = self.stack

        energy = m.energy.value
        resolution = m.resolution.value

        theta = m.get_theta(x, regular=True)

        # adding borders for smoothing:
        if not isinstance(theta, np.ndarray) or theta.ndim==0:
            theta = np.array(theta, ndmin=1)
        if theta.size==0:
            return np.ndarray(0)

        if resolution>0:
            if len(theta)>1:
                dtheta = theta[1]-theta[0]
            else:
                dtheta = resolution/5.
            blur_sigma = resolution/2.35482/dtheta

            if blur_sigma > 0.125:
                borders = int(8*blur_sigma)
                tail = np.arange(1, borders+1) * dtheta
                theta = np.hstack((theta[0] - tail[::-1],
                                   theta,
                                   theta[-1] + tail)
                                 )
        else:
            blur_sigma = 0


        periods = stack.nP # periods per group
        groupsize = stack.nGL # layers per group
        n = 1 - xi.get_optical_constants(stack._compositions,
                                         energy*1e3,
                                         stack._densities).conj()
        self.n=n

        theta_real = theta - m.offset.value

        R = _reflectivity(theta_real,
                          stack._thicknesses,
                          stack._roughnesses,
                          n,
                          12.398/energy,
                          m.polarization.value,
                          periods,
                          groupsize,
                          self.numthreads,
                          check=True)

        #R = (abs(A)**2).sum(1)

        R *= m.scale.value
        R += 10**m.background.value

        footprint = m.beam_size / np.sin(np.radians(abs(theta_real)))

        R *= np.clip(m.sample_length.value/footprint, 0.001, 1)


        if blur_sigma>0.125:
            return gaussian_filter1d(R, blur_sigma)[borders:-borders]
        else:
            return R

    def residuals(self, p, measurement, fft=True, return_sim=False):
        self.update_params(p)
        I = self.reflectivity(idm=measurement.id)
        
        if return_sim:
            return I
        
        I4 = I*measurement._q4
        # res_direct  = (I4 - R4)[m1.valid]
        res_direct  = (np.log10(I4) - np.log10(measurement.Rw))[measurement.valid]

        if fft:
            Ifft = FFT(I4[measurement.valid], measurement.fft.num)*measurement.fft.weights
            res_fourier = (Ifft - measurement.fft.value) * measurement.fft.norm
            res = np.hstack((res_direct, res_fourier.real, res_fourier.imag))
        else:
            res = res_direct
    
        mask = np.isnan(res)
        if mask.any():
            print("Warning: NaN in simulation")
            if mask.all():
                raise ValueError("Only NaN values")
            res[mask] = np.nanmean(res)*2
        return res


    def _measurements_to_DataFrame(self):
        mdata = collections.defaultdict(list)
        m_idx = []

        for measurement in self.measurements.values():
            m_idx.append(measurement.id)
            _mdata = measurement._collect_data()
            for k in _mdata:
                mdata[k].append(_mdata[k])

        return pd.DataFrame(mdata, index=m_idx)

    def _repr_html_(self):
        if not PANDAS:
            return super(Model, self).__repr__()

        title = "<h4>X-ray reflectivity model</h4>"
        html_stack = self.stack._repr_html_()
        titlem = "<h5>Measurements:</h5>"
        html_meas  = self._measurements_to_DataFrame().to_html()

        return os.linesep.join((title, html_stack, titlem, html_meas))



    def __repr__(self):
        if not PANDAS:
            return super(Stack, self).__repr__()

        str_stack = self.stack.__repr__()
        str_meas  = self._measurements_to_DataFrame().to_string()
        return os.linesep.join((str_stack, str_meas))

    def savetxt(self, path):
        with open(path, 'w') as f:
            f.write(self.__repr__())


if __name__=="__main__":
    import matplotlib.pyplot as pl
    Absorber = structure.Layer("Mo", density=8, thickness=5, roughness=0.1, name="Moly")
    Spacer = structure.Layer("B4C",  density=2.2, thickness=8., roughness=0.2, name="Boron Carbide")
    Buffer = structure.Layer("C",  3.5, 3.5, 3.5, "Carbon")

    Multilayer = structure.Group((Absorber, Spacer), 100, 1.5)

    stack = structure.Stack([Multilayer], structure.Layer("Si", 5., 3.))
    #stack = structure.Stack([Buffer], structure.Layer("Si", 5., 3.))
    stack.update()

    theta = np.linspace(0, 2, 1)
    R = np.exp(-theta)
    m = measurement.Measurement(theta, R)

    xrr = Model(stack, m)



