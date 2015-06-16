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
import numpy as np
import urllib
import sqlite3
import io
import re
import collections

#import time
electron_radius = 2.8179403e-15
avogadro = 6.022142e23
keV_A = 12.39842

DB_NAME = "materials.sqlite"
DB_PATH = os.path.join(os.path.dirname(__file__), DB_NAME)
supported_tables = ["BrennanCowan", "Chantler", "CromerLiberman", "EPDL97", 
                    "Henke", "Sasaki", "Windt"]

lsep = os.linesep


def get_components(compount, reduce_only=False):
    """
        Parse a chemical formulas of arbitrary complexity.
        
        Input:
            compount : str
                chemical formula
                e.g.: (SiO2)0.73(Na2O)0.05(K2O)0.17(CaO)0.03(Al2O3)0.02
            
            reduce_only : bool
                if True, returns sum formula only
        
        Returns:
            an elements, amount tuple where
                elements - list of element symbols
                amount   - list of quantity for each element
    """
    myfloat = lambda s: float(s) if s else 1.
    while True:
        #print compount
        parenthesis = compount.count("(")
        assert parenthesis==compount.count(")"), "Parenthesis error."
        if not parenthesis:
            break
        inner = re.findall(r'(\([^()]*\))(\d*\.*\d*)', compount)
        for (group, amount) in inner:
            fac = myfloat(amount)
            components = re.findall(r'([A-Z][a-z]?)(\d*\.*\d*)', group)
            newgroup = ["%s%g"%(k, myfloat(v)*fac) for (k,v) in components]
            compount = compount.replace(group+amount, "".join(newgroup))
    if reduce_only:
        return compount
    components = re.findall(r'([A-Z][a-z]?)(\d*\.*\d*)', compount)
    result = collections.defaultdict(float)
    for (k,v) in components:
        result[k] += myfloat(v)
    return result.keys(), result.values()

        
def store_f1f2_to_db(element, energy, f1, f2, database = DB_PATH, table = "User"):
    dbi = sqlite3.connect(database)
    cur = dbi.cursor()
    cur.execute("INSERT INTO f1f2_%s (element, energy, f1, f2) VALUES ('%s', '%f', '%f', '%f')"\
                %(table, element, energy, f1, f2))
    dbi.commit()
    dbi.close()

def get_f1f2_from_db(element, energy = None, database = DB_PATH, table = "Henke"):
    assert (table in supported_tables), \
        "input table has to be one of: %s "%supported_tables.__repr__()
    try:
        Z = int(element)
        element = get_element(Z)[2]
    except:
        Z = int(get_element(element)[-1])
    
    dbi = sqlite3.connect(database)
    cur = dbi.cursor()
    if energy==None:
        cur.execute("SELECT energy,f1,f2 FROM f1f2_" + table + " WHERE element = '%s' ORDER BY energy" %element)
        result = np.array(cur.fetchall(), dtype=float).T
        dbi.close()
        if len(result)<2:
            print("No form factors found for %s in table '%s'. Trying Henke database..." %(element, table))
            cur.execute("SELECT energy,f1,f2 FROM f1f2_Henke WHERE element = '%s' ORDER BY energy" %element)
            result = np.array(cur.fetchall(), dtype=float).T
    else:
        from scipy import interpolate
        cur.execute("SELECT energy,f1,f2 FROM f1f2_" + table + " WHERE element = '%s' ORDER BY energy" %element)
        result = np.array(cur.fetchall(), dtype=float)
        dbi.close()
        ffunc = interpolate.interp1d(result[:,0], (result[:,1], result[:,2]), bounds_error=False)
        result = ffunc(energy)
    # relativistic corrections (wrong in Sasaki, not done in Chantler)
    # X-ray data booklet (2009) Ch. 1.7 Eq 5
    if table=="Sasaki":
        result[0] += ((Z/85.455397464248506)**2.5228823203476805) 
    if table=="Chantler":
        result[0] += -((Z/82.5)**2.37)
    return result


def get_element(element, database = DB_PATH):
    dbi = sqlite3.connect(database)
    cur = dbi.cursor()
    try:
        Z = int(element)
        cur.execute("SELECT density, molar_mass, element, Z FROM elements WHERE Z = '%i'" % Z)
    except:
        cur.execute("SELECT density, molar_mass, element, Z FROM elements WHERE element = '%s'" % element)
    result=cur.fetchone()
    if not result:
        get_element_from_henke(element)
        cur.execute("SELECT density, molar_mass FROM elements WHERE element = '%s'" % element)
        result=cur.fetchone()
    rho=float(result[0])
    at_weight=float(result[1])
    symbol=str(result[2])
    Z=float(result[3])
    dbi.close()
    return rho, at_weight, symbol, Z


def store_rho_to_db(material, rho, database = DB_PATH):
    dbi = sqlite3.connect(database)
    cur = dbi.cursor()
    cur.execute("INSERT INTO densities (material, formula, density, comment) VALUES  ('%s', '%s', '%f', '%s')" % (material, material, rho, "User"))
    dbi.commit()
    dbi.close()


def get_rho_from_db(material, database = DB_PATH):
    dbi = sqlite3.connect(database)
    cur = dbi.cursor()
    if len(material)<=2 and len(get_components(material)[0])==1:
        return get_element(material)[0]
    cur.execute("SELECT density FROM densities WHERE material = '%s' OR formula = '%s' ORDER BY density_ID ASC" % (material, material))
    result=cur.fetchone()
    dbi.close()
    if not result:
        raise ValueError("Density of material '" + material + "' not found in database."+lsep)
    else:
        rho=float(result[0])
    return rho


def get_optical_constants(densities, materials, energy, database = DB_PATH, table = "Henke", feff = None):
    """
        Function to calculate the complex refractive index (delta, beta)
        of any material in the x-ray and EUV regime.

        Inputs:
            densities: list or float containing the densities of the materials
            materials: list or str (same like densities) containing the sum
                       formula of the composition of the material 
                       (e.g. 'La0.2Sr0.8MnO3')
            energy:    float or array of the X-Ray energies in eV
            database:  path to the file where the element properties and 
                       formfactors are stored.
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
    if hasattr(densities, "__iter__") and hasattr(materials, "__iter__"):
        layers=len(densities)
        if len(materials) != layers:
            raise AssertionError("densities and materials have to be lists of same length")
        delta=np.empty(layers)
        beta=np.empty(layers)
        for i in range(layers):
            elements, amount = get_components(materials[i])
            weights = np.array([get_element(element)[1]/1000. for element in elements])
            f1, f2 = np.array([get_f1f2_from_db(element, energy, database=database, table=table) for element in elements]).T
            f = np.vstack((f1,f2))
            delta[i], beta[i] = electron_radius/(2*np.pi) * (keV_A*1e-7/energy)**2 * densities[i]*1000. * avogadro * (f*amount).sum(1) / sum(amount*weights)
    elif hasattr(densities, "__float__") and type(materials)==str:
        energy = np.array(energy, ndmin=1, dtype=float)
        elements, amount = get_components(materials)
        amount = np.array(amount)
        delta = np.zeros_like(energy)
        beta = np.zeros_like(energy)
        weights = []
        for i in range(len(elements)):
            weights.append(get_element(elements[i])[1]/1000.)
            if elements[i] in feff:
                if feff[elements[i]]==0:
                    continue
                else:
                    delta += feff[elements[i]][0] * amount[i]
                    beta  += feff[elements[i]][1] * amount[i]
            else:
                f1, f2 = np.array(get_f1f2_from_db(elements[i], energy, database, table))
                delta += f1 * amount[i]
                beta  += f2 * amount[i]
        weights = np.array(weights)
        delta *= electron_radius/(2*np.pi) * (keV_A*1e-7/energy)**2 * densities*1000. * avogadro / (amount*weights).sum()
        beta  *= electron_radius/(2*np.pi) * (keV_A*1e-7/energy)**2 * densities*1000. * avogadro / (amount*weights).sum()
    else:
        raise AssertionError("densities and materials have to be lists or float and string, respectively.")
    return delta, beta


def get_attenuation_length(composition, energy, density=1, database=DB_PATH, 
                           table="Henke", feff=None):
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
            database : path to the sqlite database where the element
                       properties and formfactors are stored.
            table : set of tabulated formfactors to use.
                    can be BrennanCowan, Chantler, CromerLiberman, EPDL97, 
                           Henke, Sasaki, Windt.
                    For more details see Databases on 
                    http://www.esrf.eu/computing/scientific/dabax
                    and http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/
    """
    delta, beta = get_optical_constants(density, composition, energy, feff=feff, table=table)
    const = 10135467.657934014 # 2*eV/c/hbar
    mu = beta * const * energy
    if mu.size==1:
        mu = mu.item()
    return 1./mu