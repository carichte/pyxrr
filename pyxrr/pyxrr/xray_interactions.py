# Copyright 2013 by Carsten Richter
# Contact: carsten.richter@physik.tu-freiberg.de
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
import numpy as np
import h5py
import collections

from . import materials
from .elements import elements

#import time
electron_radius = 2.8179403e-15
avogadro = 6.022142e23
keV_A = 12.39842

DB_PATH = os.path.dirname(__file__)
F1F2_PATH = os.path.join(DB_PATH, "f1f2.h5")



_f_cache = dict()
def get_f1f2_from_db(element, energy=None, table="Henke", fallback="Henke"):
    element = elements[element]
    Z = element.number

    if element.symbol in _f_cache:
        E, f = _f_cache[element.symbol]
    else:
        with h5py.File(F1F2_PATH) as h5f:
            assert table in h5f, "table argument needs to be one of: %s"%list(h5f)
            tableh5 = h5f[table]
            if element.symbol not in tableh5:
                tableh5 = h5f[fallback]
            data = tableh5[element.symbol]
            E = data["energy"].value
            f = data["f"].value
            _f_cache[element.symbol] = E, f

    # relativistic corrections (wrong in Sasaki, not done in Chantler)
    # X-ray data booklet (2009) Ch. 1.7 Eq 5
    if table=="Sasaki":
        f += ((Z/85.455397464248506)**2.5228823203476805) 
    if table=="Chantler":
        f += -((Z/82.5)**2.37)

    if energy is not None:
        f = np.interp(energy, E, f)
        return f

    return E, f



def get_optical_constants(composition, energy, density=1, table="Henke", feff=None):
    """
        Function to calculate the complex refractive index (delta, beta)
        of any material in the x-ray and EUV regime.

        Inputs:
            composition: list or str (same like density) containing the sum
                       formula of the composition of the material 
                       (e.g. 'La0.2Sr0.8MnO3')
            energy:    float or array of the X-Ray energies in eV
            density: list or float containing the densities of the composition
            table:     set of tabulated formfactors to use.
                       can be BrennanCowan, Chantler, CromerLiberman, EPDL97, Henke, Sasaki, Windt.
                       For more details see Databases on http://www.esrf.eu/computing/scientific/dabax
                       and http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/
            feff:      a dictionary containing the effective scattering
                       amplitude of the element described by the corresponding
                       key.
    """
    if feff==None:
        feff = dict({})

    return_scalar = False
    if np.isscalar(density) and isinstance(composition, str):
        density = [density]
        composition = [composition]
        return_scalar = True

    if hasattr(density, "__iter__"):
        layers = len(density)
        if len(composition) != layers:
            raise AssertionError("length mismatch in `density` and `composition` input")
    else:
        raise AssertionError("Input for `density` or `composition` not understood.")

    n_1 = [] # refractive index - 1
    propfac = electron_radius/(2*np.pi) * (keV_A*1e-7/energy)**2 * 1000 * avogadro
    for i in range(layers):
        components, amount = materials.get_components(composition[i])
        weights = []
        f = []
        for symbol in components:
            element = elements[symbol]
            weights.append(element.atweight/1000.) # kg per mole
            if symbol in feff:
                f.append(feff[symbol])
            else:
                f.append(get_f1f2_from_db(symbol, energy, table=table))

        weights = np.array(weights)
        f = np.array(f)
        
        n_1.append(propfac * density[i] * f.T.dot(amount) / weights.dot(amount))

    if return_scalar:
        return n_1[0]

    return np.array(n_1) # delta + i*beta




def get_attenuation_length(composition, energy, density=1, table="Henke", feff=None):
    """
        Function to calculate the attenuation length 
        (inverse absorption coefficient) \mu for any material in the x-ray 
        and EUV regime.

        Inputs:
            composition : str
                containing the sum formula of the composition of the material 
                (e.g. 'La0.2Sr0.8MnO3')
            energy : float or sequence of floats
                The X-Ray energies in eV
            density : float
                density of the materials in g/cm^3
            table : set of tabulated formfactors to use.
                    can be BrennanCowan, Chantler, CromerLiberman, EPDL97, 
                           Henke, Sasaki, Windt.
                    For more details see Databases on 
                    http://www.esrf.eu/computing/scientific/dabax
                    and http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/
    """
    n_1 = get_optical_constants(density, composition, energy, feff=feff, table=table)
    const = 10135467.657934014 # 2*eV/c/hbar
    mu = n_1.imag * const * energy
    if mu.size==1:
        mu = mu.item()
    return 1./mu