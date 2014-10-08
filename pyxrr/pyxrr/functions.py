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
import xray_interactions as xi

MAX_RES=6 # maximum number of decimals in theta
lsep = os.linesep

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
        result["rho"] = xi.get_rho_from_db(result["code"])
    return result


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
                name = self.rootnames[root] + self.names["Names"][-1] \
                     + " => " + name
            else:
                name = self.rootnames[root] + name
        elif root=="grad_d":
            key = root + "_%i"%(self.dim["group"]-1)
            name = self.rootnames[root] + name
        elif self.units.has_key(root):
            name = "Meas. %i: %s (%s)"\
                   %(self.i_M-1, self.rootnames[root], self.units[root])
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
                if   fpath.lower().endswith(".njc"): data = io.load_njc(fpath)
                elif fpath.lower().endswith(".x00"): data = io.load_x00(fpath)
                elif fpath.lower().endswith(".val"): data = io.load_val(fpath)
                elif fpath.lower().endswith(".asc"): data = io.load_asc(fpath)
                elif fpath.lower().endswith(".fio"): data = io.load_fio(fpath)
                elif fpath.lower().endswith(".raw"): data = io.load_raw(fpath)
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
                    if "weighting" in props and props["weighting"] in ["z", "userdefined"] and len(data[0])>2:
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
                data[:,1]/=data[:,1].max()
                if data.shape[1]>2:
                    data[:,2]/=data[:,1].max()
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

