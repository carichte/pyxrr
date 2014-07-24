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
import os.path
import numpy as np
import urllib
import sqlite3

MAX_RES=6 # maximum number of decimals in theta
electron_radius = 2.8179403e-15
avogadro = 6.022142e23
keV_A = 12.39842

DB_NAME = "materials.sqlite"
DB_PATH = os.path.join(os.path.dirname(__file__), DB_NAME)
supported_tables = ["BrennanCowan", "Chantler", "CromerLiberman", "EPDL97", "Henke", "Sasaki", "Windt"]

lsep = os.linesep

def get_components(compount):
    complist=[]
    numlist=[]
    i=0
    while i < range(len(compount)):
	if compount[i].isupper():
	    if i==(len(compount)-1):
		complist.append(compount[i:i+1])
		numlist.append(1.)
		break
	    if compount[i+1].islower():
		complist.append(compount[i:i+2])
		if i+2==len(compount):
		    numlist.append(1.)
		    break
		x=""
		for j in range(i+2,len(compount)):
		    if compount[j].isalpha(): break
		    else: x+=compount[j]
		if x=="": x=1
		numlist.append(float(x))
		i=j
	    elif compount[i+1].isupper():
		complist.append(compount[i:i+1])
		numlist.append(1.)
		i+=1
	    elif compount[i+1].isdigit() or compount[i+1]==".":
		complist.append(compount[i:i+1])
		x=""
		for j in range(i+1,len(compount)):
		    if compount[j].isalpha(): break
		    else: x+=compount[j]
		numlist.append(float(x))
		i=j
	else: break
    return complist, numlist


def load_njc(filename):
    fd = open(filename, 'rb')
    z = fd.read(2000)
    fd.close()
    try:
        offset1 = int(z[1324])
    except:
        offset1 = int(z[1328])
    offset2 = ord(z[1348+offset1])
    offset = 1436 + offset1 + offset2
    
    fd = open(filename, 'rb')
    fd.read(offset)
    data32 = np.fromfile(file=fd, dtype=np.float32)
    fd.close()
    fd = open(filename, 'rb')
    fd.read(offset)
    data64 = np.fromfile(file=fd, dtype=np.float64)
    fd.close()
    
    start = round(data64[20], 6)
    th_range = round(data64[21], 6) - start
    step = round(data32[44], 6)
    number = int(round(th_range / step)) + 1
    time = round(data32[46], 6)
    
    angle = start + step * np.arange(number)
    intens = data32[-40-number:-40] * time
    return np.vstack((angle, intens)).T



def load_x00(filename):
    f = open(filename, 'r')
    header_length = 1
    for line in f:
        if "FirstAngle" in line:
            first = float(line.split()[1])
        if "StepWidth" in line:
            step = float(line.split()[1])
        if "TimePerStep" in line:
            time = float(line.split()[1])
        if "NrOfData" in line:
            number = int(line.split()[1])
        if "ScanData" in line:
            break
        else:
            header_length += 1
    f.close()
    
    angle = first + step * np.arange(number)
    intens = np.loadtxt(filename, skiprows=header_length) * time
    return np.vstack((angle, intens)).T



def load_fio(FILENAME):
    assert type(FILENAME)==str, AssertionError("Filename has to be of type string.")
    data=np.array([])
    fiodata=open(FILENAME, "r")
    while True:
        line = fiodata.readline()
        if not line: break
        if not line.find("!"): continue
        if "%c" in line:
            line = fiodata.readline()
            while not (line.find("!") == 0):
                if not line: break
                line = fiodata.readline()
        if "%p" in line:
            line = fiodata.readline()
            while not (line.find("!") == 0):
                if not line: break
                line = fiodata.readline()
        if "%d" in line:
            line = fiodata.readline()
            while not (line.find("!") == 0):
                if not line: break
                if "Col" in line: pass
                elif not len(data): data=np.array(line.split())
                else: data=np.vstack(( data, np.array(line.split()) ))
                line = fiodata.readline()
    fiodata.close()
    return data.astype(float)

def load_asc(FILENAME):
    fn=open(FILENAME, "r")
    zeilen=fn.readlines()
    fn.close()
    y_array=np.array([])
    for zeile in zeilen:
        if not "data" in zeile:
            if "start angle" in zeile:
                start_angle=float(zeile.split()[0])
            if "step width" in zeile:
                step_width=float(zeile.split()[0])
            if "number of steps" in zeile:
                step_num=float(zeile.split()[0])
        else:
            y_array=np.append(y_array, float(zeile.split()[0]))
    x_array=np.arange(start_angle, step_num*step_width+start_angle, step_width)
    return np.vstack((x_array, y_array)).T


def load_val(FILENAME):
    fn=open(FILENAME, "r")
    seen_intens = False
    Intensity = np.array([])
    while True:
        zeile = fn.readline()
        if not zeile: break
        if "BerPar01" in zeile:
            zeile = fn.readline()
            om_start=float(zeile)
            zeile = fn.readline()
            twth_start=float(zeile)
            zeile = fn.readline()
            twth_stop=float(zeile)
            zeile = fn.readline()
            delta_twth=float(zeile)
            zeile = fn.readline()
            zeile = fn.readline()
            zeile = fn.readline()
            zeile = fn.readline()
            steps = float(zeile)
            continue
        if "Intens" in zeile:
            seen_intens = True
            continue
        if seen_intens:
            Intensity = np.append(Intensity, float(zeile))
    fn.close()
    twotheta=np.linspace(twth_start, twth_stop, steps)
    return np.vstack((twotheta, Intensity)).T

def load_raw(FILENAME):
    # Open binary data
    fd = open(FILENAME, 'rb')
    data_32 = np.fromfile(file=fd, dtype=np.float32)
    fd.close()
    fd = open(FILENAME, 'rb')
    data_64 = np.fromfile(file=fd, dtype=np.float64)
    fd.close()
    fd = open(FILENAME, 'rb')
    data_i = np.fromfile(file=fd, dtype=np.int32)
    fd.close()
    fd = open(FILENAME, 'rb')
    data_s = np.fromfile(file=fd, dtype='a326')
    fd.close()
    # Evaluate data
    start = data_64[91]
    step = data_64[111]
    length = data_i[179]
    index = 254 + data_i[242] / 4
    twotheta = start + step * np.arange(length)
    intens = data_32[index:index+length]
    return np.vstack((twotheta, intens)).T

def read_prop_line(line, keyword):
    result = dict()
    result["sigma"] = 0
    result["periods"] = 1
    result["grad_d"] = 0
    for prop in line.replace(keyword, "").strip().split(","):
        if "=" in prop:
            thisprop=prop.partition("=")
            result[thisprop[0].strip()]=thisprop[2].strip()
    if not result.has_key("name") and result.has_key("code"):
        result["name"]=result["code"]
    if result.has_key("rho"):
        result["rho"] = float(result["rho"]) # Schichtdichten
    elif keyword in ["Ambience:", "Layer:", "Substrate:"]:
        print "Searching density of %s in database."%result["code"]
        result["rho"] = get_rho_from_db(result["code"])
    return result

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
        element = int(element)
        element = get_element(element)[2]
    except:
        pass
    dbi = sqlite3.connect(database)
    cur = dbi.cursor()
    if energy==None:
        cur.execute("SELECT energy,f1,f2 FROM f1f2_" + table + " WHERE element = '%s' ORDER BY energy" %element)
        result = np.array(cur.fetchall(), dtype=float)
        dbi.close()
        if len(result)<2:
            print("No form factors found for %s in table '%s'. Trying Henke database..." %(element, table))
            cur.execute("SELECT energy,f1,f2 FROM f1f2_Henke WHERE element = '%s' ORDER BY energy" %element)
            result = np.array(cur.fetchall(), dtype=float)
    else:
        from scipy import interpolate
        cur.execute("SELECT energy,f1,f2 FROM f1f2_" + table + " WHERE element = '%s' ORDER BY energy" %element)
        result = np.array(cur.fetchall(), dtype=float)
        dbi.close()
        ffunc = interpolate.interp1d(result[:,0], (result[:,1], result[:,2]), bounds_error=False)
        result = ffunc(energy)
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
    abbrev=str(result[2])
    Z=float(result[3])
    dbi.close()
    return rho, at_weight, abbrev, Z


def reduce_measured_data(data, new_stepping):
    """
        Deprecated
    """
    theta_curr=data[:,0].min() + new_stepping/2. #startpunkt
    columns=""
    results=None
    column_amount=len(data[0])
    zeile=np.zeros(column_amount)
    anzahl=0
    for j in range(len(data[:,0])):
        theta=data[j,0]
        if theta<(theta_curr+new_stepping/2):
            zeile[1]+=data[j,1]
            anzahl+=1
            if column_amount>2: 
                if j>0: zeile[2]=1./np.sqrt(data[j,2]**(-2) + zeile[2]**(-2)) # 1/sigma
                else: zeile[2] = data[j,2]
        else:
            zeile[0]=theta_curr
            zeile[1]/=anzahl
            if results==None: results=zeile
            else: results=np.vstack((results, zeile))
            theta_curr+=new_stepping
            anzahl=1
            zeile=data[j].copy()
    return results

def rebin_data(data, new_stepping):
    """
        Function to rebin measured data an average statistically.
        Inputs:
            data: array of measured data
                  It is assumed that:
                    first column = x
                    even columns = y
                    odd columns = error of y

            new_stepping: the new stepping for x
    """
    theta_curr=data[:,0].min() + new_stepping/2. #startpunkt
    columns=""
    results=None
    column_amount=len(data[0])
    zeile=np.zeros(column_amount)
    anzahl=0
    ind_y = (np.arange(column_amount)%2)==1
    ind_err = (np.arange(column_amount)%2)==0
    ind_err[0]=False
    for j in range(len(data[:,0])):
        theta=data[j,0]
        if theta>=(theta_curr+new_stepping/2):
            zeile[0]=theta_curr
            zeile[ind_y]/=anzahl
            zeile[ind_err]=np.sqrt(zeile[ind_err])/anzahl
            if results==None: results=zeile.copy()
            else: results=np.vstack((results, zeile))
            theta_curr+=new_stepping
            while theta>=(theta_curr+new_stepping/2.):
                zeile = np.zeros(column_amount)
                zeile[0]=theta_curr
                results=np.vstack((results, zeile))
                theta_curr+=new_stepping
            anzahl=0
            zeile = np.zeros(column_amount)
        zeile[ind_y]+=data[j,ind_y]
        zeile[ind_err]+=data[j,ind_err]**2
        anzahl+=1
    if anzahl>0:
        zeile[0]=theta_curr
        zeile[ind_y]/=anzahl
        zeile[ind_err]=np.sqrt(zeile[ind_err])/anzahl
    if results==None: results=zeile.copy()
    else: results=np.vstack((results, zeile))
    return results

def write_header(filename, header):
    """
        Write header (one line) to file
    """
    myfile = open(filename, "r")
    content = myfile.read()
    myfile.close()
    myfile = open(filename, "w")
    myfile.write(header)
    myfile.write(content)
    myfile.close()

def read_header(filename):
    """
        reads only first line of file.
    """
    myfile = open(filename, "r")
    header = myfile.readline()
    myfile.close()
    return header


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
    delta, beta = get_optical_constants(density, composition, energy)
    const = 10135467.657934014 # 2*eV/c/hbar
    mu = beta * const * energy
    if mu.size==1:
        mu = mu.item()
    return 1./mu

class ParamInput(object):
    """
        Input parameter class.
    """
    
    def __init__(self):
        self.values = dict()
        self.UniqueLayers = []
        self.names = dict(Names=[], 
                          UniqueLayers="Number of unique Layers per group")
        self.group = dict()
        self.dim = dict(rho=0, d=0, sigma=0, group=0)
        self.i_M = 0 # number of measurements
        self.rootnames = dict(rho="density (g/cm^3) ",
                              d="layer thickness (A) ",
                              sigma="interface roughness (A) ",
                              grad_d="rel. depth grading per layer (%) ",
                              energy="energy",
                              offset="angular offset",
                              background="background",
                              scale="scale",
                              resolution="resolution FWHM",
                              pol="polarization")
        
        self.units = dict(energy="keV",
                           offset="deg",
                           background="log10(I/I0)",
                           scale="1/I0",
                           resolution="deg",
                           pol="")
        self.defaults = dict(energy=8.048,
                             offset=0.,
                             background=-10.,
                             scale=1.,
                             resolution=0.,
                             pol=0.)
    
    def add(self, root, value, name=None):
        if self.dim.has_key(root):
            key = root + "_%i"%self.dim[root]
            self.dim[root] += 1
            self.group[key]  = self.dim["group"]
            if root=="rho":
                self.UniqueLayers[-1] += 1
                self.names["Names"].append(name)
            if root=="sigma":
                name = self.rootnames[root] + self.names["Names"][-1] + " => " + name
            else:
                name = self.rootnames[root] + name
        elif root=="grad_d":
            key = root + "_%i"%(self.dim["group"]-1)
            name = self.rootnames[root] + name
        elif self.units.has_key(root):
            name = "Meas. %i: %s (%s)"%(self.i_M-1, self.rootnames[root], self.units[root])
            key = root + "%i"%(self.i_M-1)
        else:
            raise ValueError("Invalid parameter name: %s"%key)
        #print key, name
        self.values[key] = value
        self.names[key]  = name
    
    def addgroup(self, periods, name=None):
        key = "N_%i"%self.dim["group"]
        self.lastN = periods
        self.values[key] = periods
        self.UniqueLayers.append(0)
        self.dim["group"] += 1
        if name!=None:
            self.names["Names"].append(name)
            self.names[key] = "Number of periods in Group %s"%name
        elif self.dim["group"] == 1:
            self.names[key] = "Number of periods in Ambience"
        else:
            self.names[key] = "Number of periods in Group %i"\
                               %(self.dim["group"]-1)
        
    
    def addmeasurement(self, **kwargs):
        self.i_M += 1
        for key in self.defaults.iterkeys():
            if kwargs.has_key(key):
                val = kwargs[key] 
            else:
                val = self.defaults[key]
            self.add(key, val)
        
    
def fetch_oc(ocdict, path, layernum):
    if not os.path.isfile(path):
        print("File not found: %s"%path)
        return
    try:
        data = np.loadtxt(path)
        Energy, delta, beta = data.T
        ocdict[layernum] = (Energy, delta, beta)
        
    except:
        print("Input format of user-given optical constants `oc' not "
              "understood for layer %i. "
              "Value: %s"%(layernum, path))
        return

def parse_parameter_file(SampleFile):
    """
        Rather complicated function to read out the parameters of a sample 
        file for pyxrr.
        
        Input: SampleFile - path to the file.
    """
    materials = []
    Parameters = ParamInput()
    total_layers = 0
    f = open(SampleFile)
    properties = f.readlines()
    f.close()
    i_M = 0
    measured_data = []
    weights = {}
    paths = []
    x_axes = []
    fit_range = {}
    lastline = ""
    oc_user = {}
    for line in properties:
        if line.strip().endswith(",") or line.strip().endswith("\\"):
            lastline+=line.strip().strip("\\")
            continue
        else:
            line = lastline + line.strip() # append continued lines
            lastline = ""
        if line.find("Ambience:") == 0:
            props=read_prop_line(line, "Ambience:")
            materials.append(props["code"])
            
            Parameters.addgroup(1)
            # Ambience thickness not relevant:
            Parameters.add("d", 0., props["name"])
            Parameters.add("grad_d", 0, props["name"])
            Parameters.add("rho", props["rho"], props["name"])
            if "oc" in props:
                fetch_oc(oc_user, props["oc"], Parameters.dim["rho"])
            
            total_layers += 1
            
        if line.find("Group:") == 0:
            props = read_prop_line(line, "Group:")
            
            if not Parameters.dim["group"]:
                raise ValueError("Define Ambience before first Group!")
            
            # Vorher Multilayer?
            if Parameters.lastN > 1:
                # Dann zuerst Rauigkeit von Letzter ML-Schicht zu erster ML-Schicht (temp)
                Parameters.add("sigma", tempsigma, tempname)
            Parameters.add("sigma", float(props["sigma"]), props["name"])
            
            
            Parameters.addgroup(int(props["periods"]), props["name"])
            Parameters.add("grad_d", float(props["grad_d"]), props["name"])
            
        if line.find("Layer:") == 0:
            props = read_prop_line(line, "Layer:")
            materials.append(props["code"])
            
            if Parameters.UniqueLayers[-1] == 0: #first layer?
                tempsigma=float(props["sigma"])
                tempname=props["name"]
            else:
                Parameters.add("sigma", float(props["sigma"]), props["name"])
            
            Parameters.add("d", float(props["d"]), props["name"])
            Parameters.add("rho", props["rho"], props["name"])
            if "oc" in props:
                fetch_oc(oc_user, props["oc"], Parameters.dim["rho"])
            total_layers+=1
        if line.find("Substrate:") == 0:
            props = read_prop_line(line, "Substrate:")
            materials.append(props["code"])
            
            # Vorher Multilayer?
            if Parameters.lastN > 1:
                # Dann zuerst Rauigkeit von Letzter ML-Schicht zu erster ML-Schicht (temp)
                Parameters.add("sigma", tempsigma, tempname)
            
            Parameters.addgroup(1) # Substrat ist kein Multilayer
            Parameters.add("sigma", float(props["sigma"]), props["name"])
            
            Parameters.add("d", 0, props["name"])
            Parameters.add("grad_d", 0, props["name"])
            Parameters.add("rho", props["rho"], props["name"])
            if "oc" in props:
                fetch_oc(oc_user, props["oc"], Parameters.dim["rho"])
            
            total_layers+=1
        if line.find("Measurement:")==0:
            props=read_prop_line(line, "Measurement:")
            print(lsep+"Parsing measurement %i..." %i_M)
            # READING MEASURED DATA FILE
            try:
                #if not os.path.isfile(props["file"]):
                #    raise IOError("File not found: %s%s"%(lsep, props["file"]))
                if os.path.isfile(props["file"]):
                    fpath = props["file"]
                    print("  Opening measured data file: %s"%fpath)
                elif os.path.isfile(os.path.join(os.path.dirname(SampleFile), props["file"])):
                    fpath = os.path.join(os.path.dirname(SampleFile), props["file"])
                    print("  Opening measured data file: %s"%fpath)
                else:
                    raise IOError("File not found: %s%s"%(lsep, props["file"]))
                if   fpath.lower().endswith(".njc"): data = load_njc(fpath)
                elif fpath.lower().endswith(".x00"): data = load_x00(fpath)
                elif fpath.lower().endswith(".val"): data = load_val(fpath)
                elif fpath.lower().endswith(".asc"): data = load_asc(fpath)
                elif fpath.lower().endswith(".fio"): data = load_fio(fpath)
                elif fpath.lower().endswith(".raw"): data = load_raw(fpath)
                else:
                    skiprows=0
                    while skiprows<200:
                        try:
                            data=np.loadtxt(fpath, skiprows=skiprows)
                            break
                        except: 
                            skiprows+=1
                    if skiprows==200:
                        raise IOError("Bad Measurement File " + fpath)
                    else:
                        print("  Ignoring first %i lines in measured data file."%(skiprows))
                paths.append(fpath)
                if props.has_key("fit_range"): 
                    limits = map(float, props["fit_range"].split("->"))
                    fit_range[i_M] = min(limits), max(limits)
                else:
                    fit_range[i_M]=0,np.inf
            except Exception as errmsg:
                if props.has_key("file"):
                    print("  Could not load measurement file %s"%props["file"])
                else: 
                    print("  Could not load measurement file")
                print("  Message: %s"%errmsg)
                data = np.array(((),())).T
                fit_range[i_M]=0,np.inf # nicht fitten
                paths.append("")
            
            # NO INF VALUES
            ind = np.isinf(data[:,1])
            if ind.any():
                print("  Infinite y-values found and discarded.")
                data = data[~ind]
            # SORT DATA
            if len(data):
                ind = data[:,0].argsort()
                data = data[ind]
            # MEASUREMENT PARAMETERS
            Parameters.addmeasurement()
            for prop_name in Parameters.defaults.iterkeys():
                if props.has_key(prop_name):
                    Parameters.add(prop_name, float(props[prop_name]))
            # WHAT KIND OF X-VALUES?
            if props.has_key("x_axis"):
                x_axes.append(props["x_axis"].lower())
            elif props.has_key("twotheta"):
                print("Warning: Keyword `twotheta' is deprecated. "
                      "Use `x_axis' instead.")
                try:
                    tth = bool(props["twotheta"])
                except:
                    raise ValueError("argument twotheta has to be boolean")
                if tth:
                    x_axes.append("twotheta")
                else:
                    x_axes.append("theta")
            else:
                x_axes.append("theta")
            
            # REBINNING
            if len(data)>0 and len(np.unique(np.diff(data[:,0]).round(MAX_RES))) > 1:
                print("  Unequal theta stepping found:")
                if props.has_key("rebin") and props["rebin"].startswith("a"):
                    print("    Data is being rebinned to maximum step width by averaging...")
                    data = rebin_data(data, np.diff(data[:,0]).max())
                else:
                    print("    Data is being rebinned to minimum step width by interpolation...")
                    from scipy.interpolate import interp1d
                    newx = np.arange(data[0,0], data[-1,0], np.diff(data[:,0]).min())
                    f_interp = interp1d(data[:,0], data[:,1])
                    if props.has_key("weighting") and props["weighting"]=="z" and len(data[0])>2:
                        f_interp2 = interp1d(data[:,0], data[:,2])
                        data = np.vstack((newx, f_interp(newx), f_interp2(newx))).T
                    else:
                        data = np.vstack((newx, f_interp(newx))).T
            # WEIGHTING
            if props.has_key("weighting"):
                weights[i_M] = props["weighting"]
            # NORMALIZATION
            if len(data)>0 and data[:,1].max()>1.:
                print("  Normalizing measured data to maximum value.")
                data[:,1]=data[:,1]/data[:,1].max()
            print("  Length of dataset %i: %i points" %(i_M, len(data)))
            print(lsep)
            i_M+=1
            measured_data.append(data)
    # NO MEASUREMENT?
    if Parameters.i_M==0:
        fit_range[i_M]=0,0 # nicht fitten
        Parameters.addmeasurement()
        for prop_name in Parameters.defaults.iterkeys():
            if props.has_key(prop_name):
                Parameters.add(prop_name, float(props[prop_name]))
        x_axes.append("theta")
        data = np.array(((),())).T # no data
        measured_data.append(data)
        
    Parameters.dim["d"] -= 1 # substrate thickness not important
    return Parameters, materials, \
           measured_data, weights, fit_range, \
           total_layers, x_axes, paths, oc_user

